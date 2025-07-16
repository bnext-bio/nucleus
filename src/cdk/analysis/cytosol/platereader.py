"""
Plate Reader Module

This module provides support for loading and analyzing data from various plate readers. Currently supported are:
+ BioTek Cytation 5
+ Revvity Envision Nexus
+ Promega Glomax Discover

Usage:
*tbd*
"""

import csv
import io
import re
import os.path
from pathlib import Path
import logging
from typing import Union, Optional, NamedTuple
from enum import Enum, auto
from ordered_set import OrderedSet

import numpy as np
import pandas as pd
from pandas.api import types as ptypes
from scipy import stats
import sklearn.metrics

# import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import timple
import timple.timedelta
import scipy.optimize
import functools
import warnings


log = logging.getLogger(__name__)

DataFile = Union[str, Path, io.StringIO]


class SteadyStateMethod(Enum):
    LOWEST_VELOCITY = (auto(),)
    MAXIMUM_VALUE = (auto(),)
    VELOCITY_INTERCEPT = auto()


class PlateReaderData(NamedTuple):
    """
    A named tuple representing plate reader data and associated plate map information.
    """

    data: pd.DataFrame
    platemap: Optional[pd.DataFrame] = None


_timple = timple.Timple()


def load_platereader_data(
    data_file: DataFile,
    platemap_file: Optional[DataFile] = None,
    platereader: Optional[str] = None,
) -> Union[PlateReaderData, pd.DataFrame]:
    """
    Load plate reader data from a file and return a DataFrame.

    This function loads platereader data from a CSV file, parsing it into a standardized format and labelling
    it with a provided plate map.

    Filenames should be formatted in a standard format: `[date]-[device]-[experiment].csv`. For
    example, `20241004-envision-dna-concentration.csv`.

    Data is loaded based on the device field in the filename, which is used to determine the appropriate reader-specific
    data parser. Currently supported readers are:
    - BioTek Cytation 5: `cytation`
    - Revvity Envision Nexus: `envision`

    Data is returned as a pandas DataFrame with the following mandatory columns:
    - `Well`: Well identifier (e.g. `A1`)
    - `Row`: Row identifier (e.g. `A`)
    - `Column`: Column identifier (e.g. `1`)
    - `Time`: Time of measurement
    - `Seconds`: Time of measurement in seconds
    - `Temperature (C)`: Temperature at time of measurement
    - `Read`: A tag describing the type of measurement (e.g. `OD600`, `Fluorescence`). The format of this field is
    currently device-specific.
    - `Data`: The measured data value

    In addition, the provided platemap will be merged to the loaded data on the `Well` column. All other columns within
    the platemap will be present in the returned dataframe.

    Args:
        data_file (str): Path to the plate reader data file.

    Returns:
        If a platemap is provided, a PlateReaderData named tuple containing the data and platemap DataFrames. Otherwise,
        just the data. If a platemap_file is provided, the returned platemap is guaranteed to be not None.

        platemap_file is not None:
            PlateReaderData: A named tuple containing the plate reader data and platemap DataFrames: (data, platemap)
        platemap_file is None:
            pd.DataFrame: DataFrame containing the plate reader data in a structured format.


    """
    if platereader is None:
        platereader = os.path.basename(data_file).lower()

    # TODO: Clean this up to use a proper platereader enum and not janky string parsing.
    if "biotek-cdk" in platereader.lower():
        data = read_biotek_cdk(data_file)
    elif "biotek" in platereader.lower():
        data = read_cytation(data_file)
    elif "cytation" in platereader.lower():
        data = read_cytation(data_file)
    elif "envision" in platereader.lower():
        data = read_envision(data_file)
    # elif filename_lower.startswith("glomax"):
    #     return read_glomax(os.path.dirname(data_file))
    else:
        raise ValueError(f"Unsupported plate reader data file: {data_file}")

    platemap = None
    if platemap_file is not None:
        platemap = read_platemap(platemap_file)
        for col in ["Row", "Column"]:
            if col in platemap:
                # Remove columns we expect might be duplicated between the platemap and
                # the data itself.
                platemap = platemap.drop(col, axis=1)

        data = data.merge(platemap, on="Well")
        return PlateReaderData(data=data, platemap=platemap)

    return data


def read_platemap(platemap_file: DataFile) -> pd.DataFrame:
    if isinstance(platemap_file, io.StringIO):
        platemap = pd.read_csv(platemap_file)
    else:
        extension = os.path.splitext(platemap_file)[1].lower()
        if extension == ".csv":
            platemap = pd.read_csv(platemap_file)
        elif extension == ".tsv":
            platemap = pd.read_table(platemap_file)
            # TODO: create test for this
        elif extension == ".xlsx":
            platemap = pd.read_excel(platemap_file)
        else:
            raise ValueError(
                f"Unsupported platemap file, use csv or xlsx: {platemap_file}"
            )

    # Remove unnamed columns from the plate map.
    platemap = platemap[
        [col for col in platemap.columns if not col.startswith("Unnamed:")]
    ]

    # Needed to make sure times are correctly converted, but we don't convert
    # floats because they get upcast to a pandas Float64Dtype() class which
    # messes up plotting.
    # platemap = platemap.convert_dtypes(convert_floating=False)

    platemap["Well"] = platemap["Well"].str.replace(
        ":", ""
    )  # Normalize well by removing : if it exists
    return platemap


# def read_glomax(data_dir: str) -> pd.DataFrame:
#     # glob over .csv files in dfpath; append to data; concatenate into one DataFrame
#     data = list()
#     for csv in glob.glob(f"{data_dir}/*.csv"):
#         df = pd.read_csv(csv)
#         df["File"] = os.path.basename(csv)
#         df["Row"] = df["WellPosition"].str.split(":").str[0]
#         df["Column"] = df["WellPosition"].str.split(":").str[1].astype(int)
#         df["Time"] = pd.to_datetime(
#             data["File"].str.replace(r".* ([0-9.]+ [0-9_]+).*", r"\1", regex=True), format="%Y.%m.%d %H_%M_%S"
#         )
#         df["WellTime"] = pd.to_timedelta(data["Timestamp(ms)"], "us")

#         data.append(df)

#     data = pd.concat(data, ignore_index=True)

#     # label different wavelengths
#     channel_map = dict(zip(data["ID"].unique(), ["A600", "Blue", "Green", "Red"]))
#     data["Channel"] = data["ID"].map(channel_map)

#     palette = dict(zip(dict.fromkeys(channel_map.values()), ["brown", "limegreen", "red", "firebrick"]))

#     # massage Time
#     data["TimeDelta"] = data["Time"] - data["Time"].min()
#     data["TimeDeltaPretty"] = data["TimeDelta"].map(
#         lambda x: "{:02d}:00".format(x.components.hours)
#     )  # {:02d}".format(x.components.hours, x.components.minutes))

#     # Get a generic data column
#     data["Data"] = data["CalculatedFlux"]
#     data["Data"].fillna(data[data["CalculatedFlux"].isna()]["OpticalDensity"], inplace=True)

#     # Label replicates
#     data["Replicate"] = data["File"].map(lambda x: re.sub(r".* OUT ([0-9]+).csv", r"\1", x))
#     data.sort_values(by=["TimeDelta", "Row", "Column"], inplace=True)

#     return data

BIOTEK_CDK_METADATA_SECTIONS = ["CDK", "Plate", "Procedure Summary"]

BIOTEK_CDK_ID_VARS = [
    "Protocol File Name",
    "Experiment File Name",
    "Plate Number",
    "Plate ID",
    "Well ID",
    "Name",
    "Well",
    "Conc/Dil type",
    "Conc/Dil",
    "Unit",
    "Time",
]

BIOTEK_CDK_PLATE_VARS = ["Plate Number", "Plate ID"]


def read_biotek_cdk(data_file: DataFile, sep="\t") -> pd.DataFrame:
    log.debug(f"Reading CDK-formatted BioTek data from {data_file}")

    with open(data_file, "r", encoding="latin1") as file:
        data_raw = file.read()

    blocks = data_raw.strip().split("\n\n")
    metadata = dict()
    dataframes = list()
    for header, section in zip(blocks[::2], blocks[1::2]):
        if header in BIOTEK_CDK_METADATA_SECTIONS:
            header_var = re.sub(r"\s+", "_", header).lower()
            log.debug(f"Found metadata: {header} -> {header_var}")

            section_metadata = [
                line.strip() for line in section.split("\n")
            ]  # Strip whitespace and break into lines
            section_metadata = [
                re.sub(r"\t+", "\t", line) for line in section_metadata
            ]  # Remove duplicated tab spaces

            section_dict = dict()
            for line in csv.reader(
                section_metadata, dialect="excel-tab", skipinitialspace=True
            ):
                field = re.sub(
                    r":$", "", line[0]
                )  # Remove the colon at the end of the field name
                section_dict[field] = line[1] if len(line) > 1 else None

            # TODO: Setting the metadata this way means we end up returning the last metadata section we saw (not per plate)
            # Return either a list, one per plate, or split it up some other way
            # Note: we are currently relying on this behavior to set the Plate ID on a data dataframe lower down.
            metadata[header_var] = section_dict

            continue

        log.debug(f"Loaded section {header}")
        data_protocol = header.split(":")[0]

        if data_protocol == "Pierce660":
            data = pd.read_table(io.StringIO(section))

            for col in data.columns:
                if col in BIOTEK_CDK_ID_VARS:
                    data[col] = data[col].ffill()

            data = data.dropna(axis=1, how="all")

            data = data.melt(
                id_vars=OrderedSet(BIOTEK_CDK_ID_VARS)
                & OrderedSet(data.columns),
                value_vars=OrderedSet(data.columns)
                - OrderedSet(BIOTEK_CDK_ID_VARS),
                var_name="Read",
                value_name="Data",
            ).reset_index()

            # data["Row"] = data["Well"].str.extract(r"[A-Z]+")
            # data["Column"] = data["Well"].str.extract(r"[0-9]+")
            data["Type"] = data["Well ID"].str.extract(r"([A-Z]+)")
            data["Sample"] = data["Well ID"].str.extract(r"([0-9]+)")

            if "Name" not in data:
                data["Name"] = data["Well ID"]
        elif data_protocol == "PURE":
            data = pd.read_table(io.StringIO(section))
            data = data.drop(columns=data.columns[1])
            data = data.dropna(axis=1, how="all")  # Remove wells with no data

            # Remove completely empty rows (usually there if the run was prematurely aborted).
            # We use data.columns[1:] because even if all wells have NaN data, the Time column will still have a time.
            # We're assuming that 'Time' is the first column.
            data = data.dropna(axis=0, subset=data.columns[1:], how="all")

            data = data.melt(
                id_vars=OrderedSet(BIOTEK_CDK_ID_VARS)
                & OrderedSet(data.columns),
                value_vars=OrderedSet(data.columns)
                - OrderedSet(BIOTEK_CDK_ID_VARS),
                var_name="Well",
                value_name="Data",
            ).reset_index(drop=True)
            log.debug(f"Data loaded with columns: {data.columns}")

            # TODO: Figure out why the plate reader adds asterisks sometimes.
            def fix_strings(x):
                if isinstance(x, str):
                    return x.replace("*", "")
                return x

            data["Data"] = data["Data"].apply(fix_strings)
            data["Data"] = pd.to_numeric(data["Data"])

            data["Time"] = pd.to_timedelta(data["Time"])

            data["Row"] = data["Well"].str.extract(r"([A-Z]+)")
            data["Column"] = data["Well"].str.extract(r"(\d+)").astype(int)
            data["Read"] = header.split(":")[1]

            if (
                "Plate Number" in metadata["cdk"]
                and metadata["cdk"]["Plate Number"] is not None
            ):
                data["Plate"] = metadata["cdk"]["Plate Number"]

            if (
                "Reading Date/Time" in metadata["cdk"]
                and metadata["cdk"]["Reading Date/Time"] is not None
            ):
                data["Clock Time"] = (
                    pd.to_datetime(metadata["cdk"]["Reading Date/Time"])
                    + data["Time"]
                )

        data.attrs["metadata"] = metadata
        dataframes.append(data)

    return pd.concat(dataframes)


def read_cytation(data_file: DataFile, sep="\t") -> pd.DataFrame:
    log.debug(f"Reading Cytation data from {data_file}")
    # read data file as long string
    data = ""
    with open(data_file, "r", encoding="latin1") as file:
        data = file.read()

    # extract indices for Proc Details, Layout
    procidx = re.search(r"Procedure Details", data)
    layoutidx = re.search(r"Layout", data)
    readidx = re.search(r"^(Read\s)?\d+(/\d+)?,\d+(/\d+)?", data, re.MULTILINE)

    # get header DataFrame
    header = data[: procidx.start()]
    header = pd.read_csv(
        io.StringIO(header), delimiter=sep, header=0, names=["key", "value"]
    )

    # get procedure DataFrame
    procedure = data[procidx.end() : layoutidx.start()]
    procedure = pd.read_csv(
        io.StringIO(procedure), skipinitialspace=True, names=range(4)
    )
    procedure = procedure.replace(np.nan, "")

    # get Cytation plate map from data_file as DataFrame
    layout = data[layoutidx.end() : readidx.start()]
    layout = pd.read_csv(io.StringIO(layout), index_col=False)
    layout = layout.set_index(layout.columns[0])
    layout.index.name = "Row"

    # iterate over data string to find individual reads
    reads = dict()

    sep = (
        r"(?:Read\s\d+:)?(?:\s\d{3}(?:/\d+)?,\d{3}(?:/\d+)?(?:\[\d\])?)?" + sep
    )

    for readidx in re.finditer(
        r"^(Read\s)?\d+(/\d+)?,\d+(/\d+)?.*\n", data, re.MULTILINE
    ):
        # for each iteration, extract string from start idx to end icx
        read = data[readidx.end() :]
        read = read[
            : re.search(
                r"(^(Read\s)?\d+,\d+|^Blank Read\s\d|Results|\Z)",
                read[1:],
                re.MULTILINE,
            ).start()
        ]
        read = pd.read_csv(
            io.StringIO(read), sep=sep, engine="python"
        ).convert_dtypes(convert_floating=False)
        reads[data[readidx.start() : readidx.end()].strip()] = read

    # create a DataFrame for each read and process, then concatenate into a large DataFrame
    # NOTE: JC 2024-05-21 - turns out, len(list(reads.items())) = 1 (one big mono table)
    read_dataframes = list()
    for name, r in reads.items():
        # filter out Cytation calculated kinetic parameters, which are cool, but don't want rn
        r = r[r.Time.str.contains(r"\d:\d{2}:\d{2}", regex=True)]

        # extract meaningful parameters from really big string
        r = r.melt(id_vars=["Time", "T°"], var_name="Well", value_name="Data")
        r["Row"] = r["Well"].str.extract(r"([A-Z]+)")
        r["Column"] = r["Well"].str.extract(r"(\d+)").astype(int)
        r["Temperature (C)"] = r["T°"]  # .str.extract(r"(\d+)").astype(float)
        r["Data"] = r["Data"].replace("OVRFLW", np.inf)
        r["Data"] = r["Data"].astype(float)
        r["Read"] = name
        r["Ex"] = r["Read"].str.extract(r"(\d+),\d+").astype(int)
        r["Em"] = r["Read"].str.extract(r"\d+,(\d+)").astype(int)
        read_dataframes.append(r)

    data = pd.concat(read_dataframes)

    # add time column to data DataFrame
    data["Time"] = pd.to_timedelta(data["Time"])
    data["Seconds"] = data["Time"].map(lambda x: x.total_seconds())

    return data[
        [
            "Well",
            "Row",
            "Column",
            "Time",
            "Seconds",
            "Temperature (C)",
            "Read",
            "Data",
        ]
    ]


def read_envision(data_file: DataFile) -> pd.DataFrame:
    # load data
    data = pd.read_csv(data_file).convert_dtypes()

    # massage Row, Column, and Well information
    data["Row"] = (
        data["Well ID"].apply(lambda s: s[0]).astype(pd.StringDtype())
    )
    data["Column"] = data["Well ID"].apply(lambda s: str(int(s[1:])))
    data["Well"] = data.apply(
        lambda well: f"{well['Row']}{well['Column']}", axis=1
    )

    data["Time"] = pd.to_timedelta(data["Time [hhh:mm:ss.sss]"])
    data["Seconds"] = data["Time"].map(lambda x: x.total_seconds())

    data["Temperature (C)"] = data["Temperature current[°C]"]

    data["Read"] = data["Operation"]

    data["Data"] = data["Result Channel 1"]

    data["Excitation (nm)"] = data["Exc WL[nm]"]
    data["Emission (nm)"] = data["Ems WL Channel 1[nm]"]
    data["Wavelength (nm)"] = (
        data["Excitation (nm)"] + "," + data["Emission (nm)"]
    )

    return data[
        [
            "Well",
            "Row",
            "Column",
            "Time",
            "Seconds",
            "Temperature (C)",
            "Read",
            "Data",
        ]
    ]


def blank_data(data: pd.DataFrame, blank_type="Blank"):
    """
    Blank data from plate reader measurements.

    Adjusts plate reader data by subtracting the value of one or more blanks at each timepoint, for each read channel.
    By default, the data will be blanked against the mean value of all wells of type "Blank".

    This function adjusts the main "Data" column in the dataframe provided, so that blanked values can be easily
    used in subsequent processing. The original (unblanked) data is available in a new 'Data_unblanked' column. The
    blank value calculated for each row of the data is present in 'Data_blank'.

    Args:
        data (pd.DataFrame): Input DataFrame containing 'Well', 'Time', 'Data', and 'Type' columns.
        blank_type (str, optional): Value in the 'Type' column to use as blank. Defaults to "Blank".

    Returns:
        pd.DataFrame: DataFrame with blanked 'Data' values and an additional 'Data_unblanked' column.

    """
    blank = (
        data[data["Type"] == blank_type]
        .groupby(["Time", "Read"])["Data"]
        .mean()
    )
    data = data.merge(
        blank, on=["Time", "Read"], suffixes=("", "_blank"), how="left"
    )

    # Check to make sure we don't have missing blanks for certain Time/Read combinations in the source data.
    # The most likely way this could happen is if the platereader "Time" isn't aligned well-to-well.
    if data["Data_blank"].isna().any():
        log.warning(
            "Not all data has a blank value; blanked data will contain NaNs."
        )

    data["Data_unblanked"] = data["Data"].copy()
    data["Data"] = data["Data"] - data["Data_blank"]

    return data


def plot_setup() -> None:
    _timple.enable()
    pd.set_option("display.float_format", "{:.2f}".format)


def _plot_timedelta(plot: sns.FacetGrid | mpl.axes.Axes) -> None:
    axes = [plot]
    if isinstance(plot, sns.FacetGrid):
        axes = plot.axes.flatten()

    for ax in axes:
        # ax.xaxis.set_major_locator(timple.timedelta.AutoTimedeltaLocator(minticks=3))
        ax.xaxis.set_major_formatter(
            timple.timedelta.TimedeltaFormatter("%h:%m")
        )
        ax.set_xlabel("Time (hours)")

    # g.set_xlabels("Time (hours)")
    # g.figure.autofmt_xdate()


def plot_plate(data: pd.DataFrame) -> sns.FacetGrid:
    g = sns.relplot(
        data=data,
        x="Time",
        y="Data",
        row="Row",
        col="Column",
        hue="Read",
        kind="line",
    )
    _plot_timedelta(g)

    g.set_ylabels("Fluorescence (RFU)")
    g.set_titles("{row_name}{col_name}")

    return g


def plot_curves_by_name(
    data: pd.DataFrame, by_experiment=True
) -> sns.FacetGrid:
    """
    Produce a basic plot of timeseries curves, coloring curves by the `Name` of the sample.

    If there are multiple different `Read`s in the data (e.g., GFP, RFP), then a subplot will be
    produced for each read. If there are multiple experiments, each experiment will be plotted separately.

    Args:
        data (pd.DataFrame): DataFrame containing plate reader data.
        by_experiment (bool, optional): If True, then each experiment will be plotted in a separate subplot.

    Returns:
        sns.FacetGrid: Seaborn FacetGrid object containing the plot.
    """
    kwargs = {}
    if "col" not in kwargs and "Experiment" in data.columns and by_experiment:
        kwargs["col"] = "Experiment"

    if "row" not in kwargs and data["Read"].unique().size > 1:
        kwargs["row"] = "Read"

    g = plot_curves(data=data, x="Time", y="Data", hue="Name", **kwargs)

    return g


def plot_curves(
    data: pd.DataFrame,
    x="Time",
    y="Data",
    hue="Name",
    labels=(None, "Fluorescence (RFU)"),
    **kwargs,
) -> sns.FacetGrid:
    """
    Plot timeseries curves from a plate reader dataset, allowing selection of the parameters to
    use for plotting and to divide the data into multiple subplots.

    This function is a thin wrapper around Seaborn `relplot`, providing sensible defaults while
    also allowing for the use of any `relplot` parameter.

    Args:
        data (pd.DataFrame): DataFrame containing plate reader data.
        x (str, optional): Column name to use for x-axis. Defaults to "Time".
        y (str, optional): Column name to use for y-axis. Defaults to "Data".
        hue (str, optional): Column name to use for color coding. Defaults to "Name".
        labels (tuple, optional): Labels for the x and y axes. Defaults to (None, "Fluorescence (RFU)").
                                  If None, use the default label (the name of the field, or a formatted time label).
        **kwargs: Additional keyword arguments passed to `sns.relplot`.

    Returns:
        sns.FacetGrid: A FacetGrid object containing the plotted data.

    """
    if "row" not in kwargs and "col" in kwargs:
        kwargs["col_wrap"] = min(data[kwargs["col"]].unique().size, 4)

    g = sns.relplot(data=data, x=x, y=y, hue=hue, kind="line", **kwargs)
    _plot_timedelta(g)

    x_label, y_label = labels
    if x_label:
        g.set_xlabels(x_label)
    if y_label:
        g.set_ylabels(y_label)

    # Set simple row and column titles, if we're faceting on row or column.
    # The join means the punctuation only gets added if we have both.
    var_len = max(
        [len(kwargs[var]) for var in ["row", "col"] if var in kwargs] + [0]
    )
    log.debug(f"{var_len=}")
    row_title = (
        f"{{row_var:>{var_len}}}: {{row_name}}" if "row" in kwargs else ""
    )
    col_title = (
        f"{{col_var:>{var_len}}}: {{col_name}}" if "col" in kwargs else ""
    )
    g.set_titles("\n".join(filter(None, [row_title, col_title])))

    return g


###
# Kinetics Analysis
# TODO: Perhaps split this out into a submodule.
###


def find_steady_state_for_well(well):
    well = well.sort_values("Time")
    pct_change = well["Data"].pct_change()
    idx_maxV = pct_change.idxmax()

    ss_idx = pct_change.loc[idx_maxV:].abs().idxmin()
    ss_time = well.loc[ss_idx, "Time"]
    ss_level = well.loc[ss_idx, "Data"]

    return pd.Series(
        {"Time_steadystate": ss_time, "Data_steadystate": ss_level}
    )


def find_steady_state(
    data: pd.DataFrame, window_size=10, threshold=0.01
) -> pd.DataFrame:
    """
    Find the steady state of the "Data" column in the provided data DataFrame.

    Args:
        data (pd.DataFrame): Input DataFrame containing 'Well', 'Time', and 'Data' columns.
        window_size (int): Size of the rolling window for calculating the rate of change.
        threshold (float): Threshold for determining steady state.

    Returns:
        pd.DataFrame: DataFrame with 'Well', 'SteadyStateTime', and 'SteadyStateLevel' columns.
    """

    result = data.groupby(["Well", "Read"]).apply(find_steady_state_for_well)
    return result


def _sigmoid(x, L, k, x0):
    return L / (1 + np.exp(-k * (x - x0)))


def _data_mean(data):
    data_mean = data.groupby("Time", as_index=False)["Data"].mean()
    data_mean["Data"] = data_mean["Data"].rolling(3, min_periods=1).mean()

    return data_mean


def kinetic_analysis_per_well(
    data: pd.DataFrame, data_column="Data"
) -> pd.DataFrame:
    data_mean = _data_mean(data)
    steadystate = find_steady_state_for_well(data_mean)

    data = data.loc[data["Time"] <= steadystate["Time_steadystate"]]
    time = data["Time"].dt.total_seconds()

    # make initial guesses for parameters
    L_initial = np.max(data[data_column])
    x0_initial = np.max(time) / 4
    k_initial = (
        np.log(L_initial * 1.1 / data[data_column] - 1) / (time - x0_initial)
    ).dropna().mean() * -1.0
    p0 = [L_initial, k_initial, x0_initial]

    # attempt fitting
    params = [0, 0, 0]
    with warnings.catch_warnings():
        warnings.simplefilter("error", scipy.optimize.OptimizeWarning)
        try:
            params, _ = scipy.optimize.curve_fit(
                _sigmoid, time, data[data_column], p0=p0
            )
        except scipy.optimize.OptimizeWarning as w:
            log.debug(f"Scipy optimize warning: {w}")
        except Exception as e:
            log.warning(f"Failed to solve: {e}")
            return None
    r_squared = sklearn.metrics.r2_score(
        data[data_column], _sigmoid(time, *params)
    )
    log.debug(f"Logistic fit R^2: {r_squared}")

    log.debug(f"{data['Well'].iloc[0]} Fitted params: {params}")

    # calculate velocities and velocity params
    v = (
        data_mean[data_column].diff()
        / data_mean["Time"].dt.total_seconds().diff()
    )
    log.debug(f"V = {v.shape}")

    maxV = v.max()
    maxV_d = data_mean.loc[v.idxmax(), data_column]
    maxV_time = data_mean.loc[v.idxmax(), "Time"]

    # calculate lag time
    lag = -maxV_d / maxV + maxV_time.total_seconds()
    lag_data = _sigmoid(lag, *params)

    # decile_upper = data[data_column].quantile(0.95)
    # decile_lower = data[data_column].quantile(0.05)

    # growth_s = (decile_upper - maxV_d) / maxV + maxV_time.total_seconds()

    # ss_time = data.loc[(data[data_column] > decile_upper).idxmax(), "Time"]
    # ss_d = data.loc[
    #     (data[data_column] > decile_upper).idxmax() :, data_column
    # ].mean()

    # kinetics = {
    #     # f"{data_column}_fit_d": y_fit,
    #     f"{data_column}_maxV": maxV,
    #     f"{data_column}_t_maxV": t_maxV,
    #     f"{data_column}_maxV_d": maxV_d,
    #     f"{data_column}_lag_s": lag,
    #     f"{data_column}_growth_s": growth_s,
    #     f"{data_column}_ss_s": ss_s,
    #     f"{data_column}_ss_d": ss_d,
    #     f"{data_column}_low_d": decile_lower,
    #     f"{data_column}_high_d": decile_upper,
    # }

    kinetics = {
        # f"{data_column}_fit_d": y_fit,
        ("Velocity", "Time"): maxV_time,
        ("Velocity", data_column): maxV_d,
        ("Velocity", "Max"): maxV,
        ("Lag", "Time"): pd.to_timedelta(lag, unit="s"),
        ("Lag", "Data"): lag_data,
        # f"{data_column}_growth_s": growth_s,
        ("Steady State", "Time"): steadystate["Time_steadystate"],
        ("Steady State", data_column): steadystate["Data_steadystate"],
        ("Fit", "L"): params[0],
        ("Fit", "k"): params[1],
        ("Fit", "x0"): params[2],
        ("Fit", "R^2"): r_squared,
    }

    return pd.Series(kinetics)
    # return kinetics


def kinetic_analysis(
    data: pd.DataFrame, group_by=["Name"], data_column="Data"
) -> pd.DataFrame:
    if data["Read"].unique().size > 1 and "Read" not in group_by:
        log.warning(
            "Kinetic analysis is not grouped on `Read`, but multiple different read types exist in the data. This is probably not what you want."
        )

    kinetics = data.groupby(group_by, sort=False).apply(
        functools.partial(kinetic_analysis_per_well, data_column=data_column)
    )
    return kinetics


def _format_timedelta(time, f="%h:%m"):
    try:
        return timple.timedelta.strftimedelta(time, f)
    except ValueError:
        return pd.NaT


def kinetic_analysis_summary(
    data: pd.DataFrame,
    kinetics: pd.DataFrame = None,
    group_by=["Name"],
    precision: float = 2,
):
    if kinetics is None:
        kinetics = kinetic_analysis(data, group_by=group_by)

    kinetics_styled = kinetics.style.format(precision=precision).format(
        _format_timedelta, subset=pd.IndexSlice[:, pd.IndexSlice[:, "Time"]]
    )
    return kinetics_styled


def plot_kinetics_by_well(
    data: pd.DataFrame,
    kinetics: pd.DataFrame,
    group_by: list[str],
    x: str = "Time",
    y: str = "Data",
    show_fit: bool = False,
    show_velocity: bool = False,
    show_mean: bool = False,
    annotate: bool = False,
    **kwargs,
):
    """
    Typical usage:

    > tk = kinetic_analysis(data=data, data_column="BackgroundSubtracted")
    > g = sns.FacetGrid(tk, col="Well", col_wrap=2, sharey=False, height=4, aspect=1.5)
    > g.map_dataframe(plot_kinetics, show_fit=True, show_velocity=True)
    """

    log.debug(f"Plotting kinetics for group: {group_by}")

    colors = sns.color_palette("Set2")

    ax = sns.scatterplot(data=data, x=x, y=y, color=colors[2], alpha=0.5)

    data_index = data.iloc[0].loc[group_by]
    kinetics = kinetics.loc[*data_index]

    log.debug(f"Data index: {data_index.values}")

    if kinetics.isna().any():
        log.info(f"Kinetics information not available for {data_index}.")
        return

    if show_mean:
        log.debug("Plotting data mean")
        sns.scatterplot(
            data=_data_mean(data), x=x, y=y, color=colors[4], alpha=0.5, ax=ax
        )

    # ax_ylim = (
    #     ax.get_ylim()
    # )  # Use this to run lines to bounds later, then restore them before returning.

    if show_fit:
        L = kinetics["Fit", "L"]
        k = kinetics["Fit", "k"]
        x0 = kinetics["Fit", "x0"]
        sns.lineplot(
            x=data["Time"],
            y=_sigmoid(data["Time"].dt.total_seconds(), L, k, x0),
            linestyle="--",
            color=colors[3],
            # alpha=0.5,
            ax=ax,
        )
    #     sns.lineplot(data=data, x=x, y=y, linestyle="--", c="red", alpha=0.5)

    # Max Velocity
    # maxV_x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    maxV_y = (
        kinetics["Velocity", "Max"]
        * (data["Time"] - kinetics["Velocity", "Time"]).dt.total_seconds()
        + kinetics["Velocity", "Data"]
    )

    sns.lineplot(
        x=data["Time"].loc[(maxV_y > 0) & (maxV_y < data[y].max())],
        y=maxV_y[(maxV_y > 0) & (maxV_y < data[y].max())],
        linestyle="--",
        color=colors[1],
        ax=ax,
    )

    maxV = kinetics["Velocity", "Max"]
    maxV_s = kinetics["Velocity", "Time"]
    maxV_d = kinetics["Velocity", "Data"]

    # decile_upper = summary[f"{y}_high_d"]
    # decile_lower = summary[f"{y}_low_d"]
    # ax.vlines(
    #     lag,
    #     ymin=ax_ylim[0],
    #     ymax=decile_lower,
    #     colors=colors[2],
    #     linestyle="--",
    # )

    # Time to Steady State
    ss_s = kinetics["Steady State", "Time"]
    ax.axvline(ss_s, c=colors[3], linestyle="--")

    # # Range
    # ax.axhline(decile_upper, c=colors[7], linestyle="--")
    # ax.axhline(decile_lower, c=colors[7], linestyle="--")

    if annotate:
        # Plot the text annotations on the chart
        ax.annotate(
            f"$V_{{max}} =$ {maxV:.2f} u/s",
            (maxV_s, maxV_d),
            xytext=(24, 0),
            textcoords="offset points",
            arrowprops={"arrowstyle": "->"},
            ha="left",
            va="center",
            c="black",
        )

        f = timple.timedelta.TimedeltaFormatter("%h:%m")
        lag_label = f.format_data(
            timple.timedelta.timedelta2num(kinetics["Lag", "Time"])
        )
        ax.annotate(
            f"$t_{{lag}} =$ {lag_label}",
            (kinetics["Lag", "Time"], kinetics["Lag", "Data"]),
            xytext=(12, 0),
            textcoords="offset points",
            ha="left",
            va="center",
        )

        ss_label = f.format_data(
            timple.timedelta.timedelta2num(kinetics["Steady State", "Time"])
        )
        ax.annotate(
            f"$t_{{steady state}} =$ {ss_label}",
            (
                kinetics["Steady State", "Time"],
                kinetics["Steady State", "Data"],
            ),
            xytext=(0, -12),
            textcoords="offset points",
            ha="left",
            va="top",
        )

    # Velocity
    if show_velocity:
        # TODO: This is currently broken due to rolling calculation and its effect on bounds.
        # Show a velocity sparkline over the plot
        velocity = (
            data.transform({y: "diff", x: lambda x: x}).rolling(5).mean()
        )
        velocity[y] = velocity[y]
        # velocity_ax = ax.secondary_yaxis(location="right",
        #                                  functions=(lambda x: pd.Series(x).rolling(5).mean().values, lambda x: x))
        velocity_ax = ax.twinx()
        sns.lineplot(data=velocity, x=x, y=y, alpha=0.5, ax=velocity_ax)
        velocity_ax.set_ylabel("$V (u/s)$")
        velocity_ax.set_ylim((0, velocity[y].max() * 2))

    # ax.set_ylim(ax_ylim)

    _plot_timedelta(ax)


def plot_kinetics(
    data: pd.DataFrame,
    group_by=["Name", "Read"],
    kinetics: pd.DataFrame = None,
    **kwargs,
):
    if kinetics is None:
        kinetics = kinetic_analysis(data, group_by=group_by)

    if "col_wrap" not in kwargs:
        kwargs["col_wrap"] = 2

    g = sns.FacetGrid(data, col=group_by[0], height=4, aspect=1.5, **kwargs)
    g.map_dataframe(
        plot_kinetics_by_well,
        kinetics=kinetics,
        group_by=group_by,
        show_fit=True,
        show_velocity=False,
        show_mean=True,
        annotate=True,
    )
    g.set_ylabels("Fluorescence (RFU)")


def plot_steadystate(data: pd.DataFrame, x="Name", **kwargs):
    steady_state = find_steady_state(data).reset_index()
    data_with_steady_state = steady_state.merge(
        data, on=["Well", "Read"], how="left"
    )

    # return data_with_steady_state

    if "col_wrap" not in kwargs and "col" in kwargs:
        kwargs["col_wrap"] = 2

    g = sns.catplot(
        data=data_with_steady_state,
        x=x,
        y="Data_steadystate",
        kind="bar",
        height=4,
        aspect=1.5,
        sharex=False,
        **kwargs,
    )

    g.set_xticklabels(rotation=90)
    g.set_ylabels("Steady State Fluorescence (RFU)")
    return g


def compute_standard_curve(
    data: pd.DataFrame, sc_type="STD", include_mean=True
):
    std = data[data["Type"] == sc_type]

    def curve(group):
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            group["Conc/Dil"], group["Data"]
        )
        return pd.Series(
            {
                "slope": slope,
                "intercept": intercept,
                "R": r_value,
                "R^2": r_value**2,
            }
        )

    curves = (
        std.groupby(["Type", "Read"])
        .apply(curve, include_groups=False)
        .reset_index()
    )

    if include_mean:
        mean = (
            std.groupby(["Type"])
            .apply(curve, include_groups=False)
            .reset_index()
        )
        mean["Read"] = "Mean"
        curves = pd.concat([curves, mean])

    if "Unit" in curves:
        curves["Unit"] = std["Unit"].iloc[0]

    return curves


def compute_concentration(
    data: pd.DataFrame,
    sc_type="STD",
    sample_type="SPL",
    dilution_col="Conc/Dil",
    unit_col="Unit",
):
    curves = compute_standard_curve(data, sc_type)

    conc = pd.merge(
        data[data["Type"] == sample_type],
        curves,
        on=["Read"],
        how="left",
        suffixes=("", "_y"),
    )
    conc["Concentration"] = (conc["Data"] - conc["intercept"]) / conc["slope"]

    if dilution_col in data.columns:
        conc["Original Concentration"] = (
            conc["Concentration"] * conc[dilution_col]
        )

    return conc


def plot_standard_curve(data: pd.DataFrame, sc_type="STD", **kwargs):
    curves = compute_standard_curve(data, sc_type)

    x_min = data.loc[data["Type"] == sc_type, "Conc/Dil"].min()
    x_max = data.loc[data["Type"] == sc_type, "Conc/Dil"].max()

    lines = list()
    for i, curve in curves.iterrows():
        x = np.linspace(x_min, x_max)
        y = x * curve["slope"] + curve["intercept"]
        df = pd.DataFrame(dict(x=x, y=y))
        for c in curve.index:
            df[c] = curve[c]
        lines.append(df)

    ax = sns.lineplot(
        data=pd.concat(lines),
        x="x",
        y="y",
        hue="Read",
        linestyle="--",
        **kwargs,
    )

    if ax not in kwargs:
        kwargs["ax"] = ax

    kwargs["legend"] = False
    sns.scatterplot(
        data=data[data["Type"] == sc_type],
        x="Conc/Dil",
        y="Data",
        hue="Read",
        **kwargs,
    )

    unit = ""
    if "Unit" in data:
        unit = data["Unit"].iloc[0]
        unit = f"({unit})"
    ax.set_xlabel(f"Concentration {unit}")

    ax.set_ylabel("Relative Fluorescence Units (RFU)")

    return ax


def plot_concentration(
    data: pd.DataFrame, x="Name", y=None, hue="Read", **kwargs
):
    conc = compute_concentration(data)

    y_column = y if y is not None else "Concentration"
    if "Original Concentration" in conc.columns:
        y_column = "Original Concentration"

    ax = sns.stripplot(data=conc, x=x, y=y_column, hue=hue, **kwargs)

    means = conc.groupby(["Name"])["Concentration"].mean().reset_index()

    sns.boxplot(
        data=means,
        x="Name",
        y="Concentration",
        showmeans=True,
        meanline=True,
        meanprops={"color": "k", "ls": "--", "lw": 1, "alpha": 0.5},
        medianprops={"visible": False},
        whiskerprops={"visible": False},
        zorder=10,
        showfliers=False,
        showbox=False,
        showcaps=False,
        fill=False,
        legend=True,
        dodge=False,
        # ax=ax
    )

    # if "Unit" in conc:
    #     ax.set_ylabel(f"Concentration ({conc['Unit'].iloc[0]})")

    return ax


def export_data(data: pd.DataFrame, output_file: Path):
    """
    Exports platereader data to a CSV.

    Args:
        data (pd.DataFrame): The plate reader data to export.
        output_file (os.path.Path): The output file.

    Returns:
        None
    """
    data.to_csv(output_file, index=False)


def export_kinetics(
    kinetics: pd.DataFrame, output_file: Path, platemap: pd.DataFrame = None
):
    """
    Exports the kinetics analysis dataframe to a CSV file.

    Args:
        kinetics (pd.DataFrame): Kinetic results (the output of `kinetics_analysis()`).
        output_file (os.path.Path): Path to the output CSV file.
        platemap (pd.Dataframe): Platemap. If this is provided, it will be merged into the kinetics results to provide full labels.

    Returns:
        None
    """

    df = kinetics.copy()
    kinetics_index = kinetics.index.names
    df.columns = df.columns.map("_".join)
    df = df.reset_index()

    for i, col in enumerate(df.columns):
        if ptypes.is_timedelta64_dtype(df[col]):
            df[col] = df[col].dt.total_seconds()
            df = df.rename(columns={col: f"{col} (s)"})

    if platemap is not None:
        df = df.merge(platemap, how="left", on=kinetics_index)

    df.to_csv(output_file, index=False)


def merge_plates(data: pd.DataFrame, plates: list[str] = None) -> pd.DataFrame:
    """Merge multiple plates in a timeseries into one plate, adjusting times.

    Where there is more than one plate in a dataframe (specified by the `Plate` column),
    merge them together into one continous timeseries. This is useful if, for example, the
    plate reader was stopped and restarted, so the two plates are really one.

    `Time` offsets after the first plate are adjusted based on read start time of each plate.

    To work, the `data` dataframe needs several columns:
    + `Plate`
    + `Clock Time`

    Args:
        data (pd.DataFrame): the plate reader data
        plates (list[str]): a list of the plates to merge. By default, all plates will be merged in the order they appear.
    """

    if plates is None:
        plates = data["Plate"].unique()

    if len(plates) <= 1:
        log.warning(f"Data does not have multiple plates to merge: {plates}")
        return data

    plate_data = [data[data["Plate"] == plates[0]]]

    start_time = data.loc[data["Plate"].isin(plates), "Clock Time"].min()
    log.debug(f"Start time: {start_time}")

    for plate in plates[1:]:
        p = data[data["Plate"] == plate].copy()
        p["Time"] = p["Clock Time"] - start_time
        plate_data.append(p)

    merged_data = pd.concat(plate_data)
    merged_data["Plate"] = plates[0]

    return merged_data
