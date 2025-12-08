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
import math
from itertools import combinations
import re

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import timple
import timple.timedelta
import scipy.optimize
import warnings

TIME_COLUMN_NAME = "Time"
DATA_COLUMN_NAME = "Data"  # fluorescence data
READ_COLUMN_NAME = "Read"
DEFAULT_FIT_FUNCTION_NAME = "sigmoid_drift"
STEADY_STATE = "Steady State"

# Track which column downstream functions should use by default
ACTIVE_DATA_COLUMN_NAME = DATA_COLUMN_NAME

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

def set_active_data_column(col: str):
    global ACTIVE_DATA_COLUMN_NAME
    ACTIVE_DATA_COLUMN_NAME = col


def get_active_data_column() -> str:
    return ACTIVE_DATA_COLUMN_NAME

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
    TIME_COLUMN_NAME,
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
        data = None

        if data_protocol == "Pierce660" or data_protocol == "Endpoint":
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
                var_name=READ_COLUMN_NAME,
                value_name=DATA_COLUMN_NAME,
            ).reset_index()

            # data["Row"] = data["Well"].str.extract(r"[A-Z]+")
            # data["Column"] = data["Well"].str.extract(r"[0-9]+")
            if "Well ID" in data:
                data["Type"] = data["Well ID"].str.extract(r"([A-Z]+)")
                data["Sample"] = data["Well ID"].str.extract(r"([0-9]+)")

            if "Name" not in data:
                if "Well ID" in data:
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
                value_name=DATA_COLUMN_NAME,
            ).reset_index(drop=True)
            log.debug(f"Data loaded with columns: {data.columns}")

            # TODO: Figure out why the plate reader adds asterisks sometimes.
            def fix_strings(x):
                if isinstance(x, str):
                    if x == "OVRFLW":
                        return np.nan
                    return x.replace("*", "")
                return x

            data[DATA_COLUMN_NAME] = data[DATA_COLUMN_NAME].apply(fix_strings)
            data[DATA_COLUMN_NAME] = pd.to_numeric(data[DATA_COLUMN_NAME])

            data[TIME_COLUMN_NAME] = pd.to_timedelta(data[TIME_COLUMN_NAME])

            data["Row"] = data["Well"].str.extract(r"([A-Z]+)")
            data["Column"] = data["Well"].str.extract(r"(\d+)").astype(int)
            data[READ_COLUMN_NAME] = header.split(":")[1]

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
                    + data[TIME_COLUMN_NAME]
                )

        if data is None:
            log.warning(
                f"No data loaded: section potentially unrecognized: {data_protocol}"
            )
            continue

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
        r = r.melt(
            id_vars=[TIME_COLUMN_NAME, "T°"],
            var_name="Well",
            value_name=DATA_COLUMN_NAME,
        )
        r["Row"] = r["Well"].str.extract(r"([A-Z]+)")
        r["Column"] = r["Well"].str.extract(r"(\d+)").astype(int)
        r["Temperature (C)"] = r["T°"]  # .str.extract(r"(\d+)").astype(float)
        r[DATA_COLUMN_NAME] = r[DATA_COLUMN_NAME].replace("OVRFLW", np.inf)
        r[DATA_COLUMN_NAME] = r[DATA_COLUMN_NAME].astype(float)
        r[READ_COLUMN_NAME] = name
        r["Ex"] = r[READ_COLUMN_NAME].str.extract(r"(\d+),\d+").astype(int)
        r["Em"] = r[READ_COLUMN_NAME].str.extract(r"\d+,(\d+)").astype(int)
        read_dataframes.append(r)

    data = pd.concat(read_dataframes)

    # add time column to data DataFrame
    data[TIME_COLUMN_NAME] = pd.to_timedelta(data[TIME_COLUMN_NAME])
    data["Seconds"] = data[TIME_COLUMN_NAME].map(lambda x: x.total_seconds())

    return data[
        [
            "Well",
            "Row",
            "Column",
            TIME_COLUMN_NAME,
            "Seconds",
            "Temperature (C)",
            READ_COLUMN_NAME,
            DATA_COLUMN_NAME,
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

    data[TIME_COLUMN_NAME] = pd.to_timedelta(data["Time [hhh:mm:ss.sss]"])
    data["Seconds"] = data[TIME_COLUMN_NAME].map(lambda x: x.total_seconds())

    data["Temperature (C)"] = data["Temperature current[°C]"]

    data[READ_COLUMN_NAME] = data["Operation"]

    data[DATA_COLUMN_NAME] = data["Result Channel 1"]

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
            TIME_COLUMN_NAME,
            "Seconds",
            "Temperature (C)",
            READ_COLUMN_NAME,
            DATA_COLUMN_NAME,
        ]
    ]


def blank_data(data: pd.DataFrame, blank_type="Blank"):
    """
    Blank data from plate reader measurements.

    Adjusts plate reader data by subtracting the value of one or more blanks at each timepoint, for each read channel.
    By default, the data will be blanked against the mean value of all wells of type "Blank".

    This function adjusts the main DATA_COLUMN_NAME column in the dataframe provided, so that blanked values can be easily
    used in subsequent processing. The original (unblanked) data is available in a new 'Data_unblanked' column. The
    blank value calculated for each row of the data is present in 'Data_blank'.

    Args:
        data (pd.DataFrame): Input DataFrame containing 'Well', 'Time', 'Data', and 'Type' columns.
        blank_type (str, optional): Value in the 'Type' column to use as blank. Defaults to "Blank".

    Returns:
        pd.DataFrame: DataFrame with blanked 'Data' values and an additional 'Data_unblanked' column.

    """
    data_column = get_active_data_column()
    blank = (
        data[data["Type"] == blank_type]
        .groupby([TIME_COLUMN_NAME, READ_COLUMN_NAME])[data_column]
        .mean()
    )
    data = data.merge(
        blank,
        on=[TIME_COLUMN_NAME, READ_COLUMN_NAME],
        suffixes=("", "_blank"),
        how="left",
    )

    # Check to make sure we don't have missing blanks for certain Time/Read combinations in the source data.
    # The most likely way this could happen is if the platereader TIME_COLUMN_NAME isn't aligned well-to-well.
    if data["Data_blank"].isna().any():
        log.warning(
            "Not all data has a blank value; blanked data will contain NaNs."
        )

    data["Data_unblanked"] = data[data_column].copy()
    data[data_column] = data[data_column] - data["Data_blank"]

    return data


def plot_setup() -> None:  # TODO:make this not have to be called explicitly...
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
        x=TIME_COLUMN_NAME,
        y=get_active_data_column(),
        row="Row",
        col="Column",
        hue=READ_COLUMN_NAME,
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

    if "row" not in kwargs and data[READ_COLUMN_NAME].unique().size > 1:
        kwargs["row"] = READ_COLUMN_NAME

    g = plot_curves(
        data=data, x=TIME_COLUMN_NAME, y=get_active_data_column(), hue="Name", **kwargs
    )

    return g


def plot_curves(
    data: pd.DataFrame,
    x=TIME_COLUMN_NAME,
    y=None,
    hue="Name",
    labels=(None, "Fluorescence (RFU)"),
    show_plot=False,
    **kwargs,
) -> sns.FacetGrid:
    """
    Plot timeseries curves from a plate reader dataset, allowing selection of the parameters to
    use for plotting and to divide the data into multiple subplots.

    This function is a thin wrapper around Seaborn `relplot`, providing sensible defaults while
    also allowing for the use of any `relplot` parameter.

    Args:
        data (pd.DataFrame): DataFrame containing plate reader data.
        x (str, optional): Column name to use for x-axis. Defaults to TIME_COLUMN_NAME.
        y (str, optional): Column name to use for y-axis. Defaults to ACTIVE_DATA_COLUMN_NAME.
        hue (str, optional): Column name to use for color coding. Defaults to "Name".
        labels (tuple, optional): Labels for the x and y axes. Defaults to (None, "Fluorescence (RFU)").
                                  If None, use the default label (the name of the field, or a formatted time label).
        **kwargs: Additional keyword arguments passed to `sns.relplot`.

    Returns:
        sns.FacetGrid: A FacetGrid object containing the plotted data.

    """
    if y is None:
        y = get_active_data_column()

    col = kwargs.get("col", None)
    row = kwargs.get("row", None)
    if (
        data[READ_COLUMN_NAME].nunique() > 1
            and col != READ_COLUMN_NAME
            and row != READ_COLUMN_NAME
    ):
        log.warning(
            "Multiple different read types exist in the data. But plotting is not grouped by \"Read\". "
            "Defaulted to plotting all Read types \n"
            "In the future, either explicitly add col=\"Read\", or provide a dataset with only one read type \n"
            "You can set one read type by setting data= data[data.Read == \"<YOUR DESIRED GAIN>\"], "
            f"where <> is replaced by one of the following:  {data.READ_COLUMN_NAME.unique()} "
        )
        kwargs["col"] = READ_COLUMN_NAME

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

    if show_plot:
        plt.show()
    return g


###
# Kinetics Analysis
# TODO: Perhaps split this out into a submodule.
###


def estimate_steady_state_for_well(well, window=3, data_col=None):
    """
    Estimates steady state by finding where percent change stabilizes.
    :param well:
    :param window:
    :param data_col:
    :return:
    """
    if data_col is None:
        data_col = get_active_data_column()

    well = well.sort_values(TIME_COLUMN_NAME)
    well_mean = _data_mean(well, window=window, data_col=data_col)
    pct_change = well_mean[data_col].pct_change()
    idx_maxV = pct_change.idxmax()

    ss_idx = pct_change.loc[idx_maxV:].abs().idxmin()
    ss_time = well_mean.iloc[ss_idx][TIME_COLUMN_NAME]
    ss_level = well_mean.iloc[ss_idx][data_col]

    return pd.Series(
        {
            f"{TIME_COLUMN_NAME}_steadystate": ss_time,
            f"{data_col}_steadystate": ss_level,
        }
    )


def find_steady_state(
    data: pd.DataFrame,
    group_by=["Name", READ_COLUMN_NAME],
    window=3,
    data_column=None,
    fit_function_name=DEFAULT_FIT_FUNCTION_NAME,
) -> pd.DataFrame:
    """
    Find the steady state of the ACTIVE_DATA_COLUMN_NAME column in the provided data DataFrame.

    Args:
        data (pd.DataFrame): Input DataFrame containing 'Well', 'Time', and 'Data' columns.
        window_size (int): Size of the rolling window for calculating the rate of change.
        threshold (float): Threshold for determining steady state.

    Returns:
        pd.DataFrame: DataFrame with 'Well', 'SteadyStateTime', and 'SteadyStateLevel' columns.
    """
    if data_column is None:
        data_column = get_active_data_column()

    if (
        data[READ_COLUMN_NAME].unique().size > 1
        and READ_COLUMN_NAME not in group_by
    ):
        log.warning(
            "Data is not grouped on `Read`, but multiple different read types exist in the data. This is probably not what you want."
        )

    if fit_function_name == "sigmoid":
        result = data.groupby(group_by, sort=False).apply(
            estimate_steady_state_for_well, data_col=data_column
        )
    else:
        kinetics_fits = kinetic_analysis(
            data, group_by=group_by, fit_function_name=fit_function_name, data_column=data_column
        )
        result = kinetics_fits.loc[
            :, [(STEADY_STATE, TIME_COLUMN_NAME), (STEADY_STATE, data_column)]
        ]
        # Keep labels the same as the old method ("Time_steadystate" and "Data_steadystate")
        # TODO: eventually make this all be the same column names as from kinetics
        result[f"{TIME_COLUMN_NAME}_steadystate"] = kinetics_fits.loc[
            :, (STEADY_STATE, TIME_COLUMN_NAME)
        ]
        result[f"{data_column}_steadystate"] = kinetics_fits.loc[
            :, (STEADY_STATE, data_column)
        ]
    return result


def _sigmoid(x, L, k, x0):
    return L / (1 + np.exp(-k * (x - x0)))


def _sigmoid_drift(t, L, k, t0, b, tau):
    """

    :param t:
    :param L:
    :param k:
    :param t0:
    :param b: drift rate
    :param tau: drift at onset
    :return:
    """
    return _sigmoid(t, L, k, t0) + b * (t - tau)

def l2_loss_sigmoid(params, t, y_data, use_drift = DEFAULT_FIT_FUNCTION_NAME):
    if use_drift == "sigmoid_drift":
        # params: [L, k, t0, b, tau]
        L, k, t0, b, tau = params
        y_model = _sigmoid_drift(t, L, k, t0, b, tau)
    else:
        # params: [L, k, x0]
        L, k, t0 = params
        y_model = _sigmoid(t, L, k, t0)
    # L2 (least squares) loss
    return np.sum((y_data - y_model) ** 2)

def  logistic_time_to_fraction(A, k, t0, alpha=0.95):
    """
    Compute the time when a logistic curve reaches a fraction alpha of its asymptote A.

    Logistic model:
        y(t) = A / (1 + exp(-k * (t - t0)))

    Parameters
    ----------
    A : float
        Asymptote (maximum value).
    k : float
        Growth rate.
    t0 : float
        Inflection point (time at which y = A/2).
    alpha : float
        Fraction of A to solve for (0 < alpha < 1).
        Example: alpha=0.95 gives time when y(t) = 0.95*A.

    Returns
    -------
    t_alpha : float
        Time when y(t) = alpha * A.
    """
    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1 (exclusive).")
    if k == 0:
        return 0
    return t0 - (1.0 / k) * np.log((1.0 / alpha) - 1.0)


FIT_FUNCTIONS = {"sigmoid": _sigmoid, "sigmoid_drift": _sigmoid_drift}


def _data_mean(data, window=3, data_col=None) -> pd.DataFrame:
    if data_col is None:
        data_col = get_active_data_column()

    data_mean = data.groupby(TIME_COLUMN_NAME, as_index=False)[data_col].mean()
    data_mean[data_col] = (
        data_mean[data_col].rolling(window, min_periods=1).mean()
    )

    return data_mean


def kinetic_analysis_per_well(
    data: pd.DataFrame,
    data_column=None,
    fit_function_name=DEFAULT_FIT_FUNCTION_NAME,
    group_keys=["Name", "Read", "Well"],
) -> pd.Series:
    """
    Perform kinetic analysis on a single well or group of data.

    Uses either sigmoid or sigmoid_drift fits.

    Args:
        data: DataFrame containing time series for a single well/group.
        group_keys: Values of the groupby keys for logging.
        data_column: Column to fit/analyze.
        fit_function_name: 'sigmoid' or 'sigmoid_drift'.

    Returns:
        pd.Series: Multi-index series with kinetics results.
    """
    if data_column is None:
        data_column = get_active_data_column()

    if data[data_column].dropna().empty:
        log.warning(
            f"{group_keys} had no data for {data_column}. Possibly from overflow errors. No fit could be made."
        )
        # Return a NaN-filled series to keep table structure consistent
        return pd.Series(
            {
                ("Velocity", TIME_COLUMN_NAME): np.nan,
                ("Velocity", data_column): np.nan,
                ("Velocity", "Max"): np.nan,
                ("Lag", TIME_COLUMN_NAME): np.nan,
                ("Lag", data_column): np.nan,
                (STEADY_STATE, TIME_COLUMN_NAME): np.nan,
                (STEADY_STATE, data_column): np.nan,
                ("Fit", "params"): [],
                ("Fit", "R^2"): np.nan,
                ("Fit", "drift"): np.nan,
                ('Fit', 'R^2'): np.nan,
                ("Fit", "drift"): np.nan,
                ("Fit", 'good_fit'): False
            }
        )
    fit_function = FIT_FUNCTIONS[fit_function_name]

    steadystate = estimate_steady_state_for_well(data, data_col=data_column)

    if fit_function_name == "sigmoid":
        data = data.loc[
            data[TIME_COLUMN_NAME]
            <= steadystate[f"{TIME_COLUMN_NAME}_steadystate"]
        ]
    time = data[TIME_COLUMN_NAME].dt.total_seconds()/3600.0

    data_mean = _data_mean(data, data_col=data_column)

    # make initial guesses for parameters
    L_initial = np.max(data[data_column])

    # the range of conditions we see in the data.
    data_mean_v = data_mean[data_column].rolling(3).mean().diff()
    timestep = (
        data_mean[TIME_COLUMN_NAME]
        .dt.total_seconds()
        .sort_values()
        .diff()
        .mean()/3600.0
    )

    k_initial = data_mean_v.max() / timestep / L_initial * 4
    x0_initial = (
        data_mean[TIME_COLUMN_NAME]
        .dt.total_seconds()
        .loc[data_mean_v.idxmax()]/3600.0
    )

    p0 = [L_initial, k_initial, x0_initial]
    log.debug(
        f"Estimated initial parameters: L={L_initial:.4f}, k={k_initial:.4f}, x0={x0_initial:.4f}"
    )

    # attempt fitting
    params = [0, 0, 0]
    bounds = [(0, L_initial*10), (0, 10000), (0,  np.max(time))] #,
    if fit_function_name == "sigmoid_drift":
        params = params + [0, 0]
        p0 = p0 + [k_initial / 2, 4]
        # params: [L, k, t0, b, tau]
        bounds = bounds + [(-np.inf, np.inf), (0, np.max(time))]

    with warnings.catch_warnings():
        warnings.simplefilter("error", scipy.optimize.OptimizeWarning)
        try:
            # params, _ = scipy.optimize.curve_fit(
            #     fit_function, time, data[data_column], p0=p0
            # )
            fit_results = scipy.optimize.minimize(l2_loss_sigmoid, x0=p0,
                                             args=(time, data[data_column],
                                                   fit_function_name),
                                             method='L-BFGS-B',
                                             bounds=bounds)
            params = fit_results.x
        except scipy.optimize.OptimizeWarning as w:
            log.debug(f"Scipy optimize warning: {w}")
        except Exception as e:
            log.warning(f"Failed to solve for group {group_keys}: {e}")
            return None
    r_squared = sklearn.metrics.r2_score(
        data[data_column], fit_function(time, *params)
    )
    log.debug(f"Logistic fit R^2: {r_squared}")

    log.debug(f"{group_keys} Fitted params: {params}")

    # calculate velocities and velocity params
    v = (
        data_mean[data_column].diff()
        / data_mean[TIME_COLUMN_NAME].dt.total_seconds().diff()/3600.0
    )
    log.debug(f"V = {v.shape}")

    maxV = params[1] * params[0]/ 4 #v.max()
    maxV_d = params[0] / 2 #data_mean.loc[v.idxmax(), data_column]
    maxV_time = pd.to_timedelta(params[2], unit="h") #data_mean.loc[v.idxmax(), TIME_COLUMN_NAME]

    # calculate lag time
    lag = -maxV_d / maxV + maxV_time.total_seconds()/3600.0
    lag_data = fit_function(lag, *params)

    kinetics = {
        # f"{data_column}_fit_d": y_fit,
        ("Velocity", TIME_COLUMN_NAME): maxV_time,
        ("Velocity", data_column): maxV_d,
        ("Velocity", "Max"): maxV,  # max slope, but not k which is just steepness. #maxV,
        ("Lag", TIME_COLUMN_NAME): pd.to_timedelta(lag, unit="h"),
        ("Lag", data_column): lag_data,
        # f"{data_column}_growth_s": growth_s,
        (STEADY_STATE, TIME_COLUMN_NAME): steadystate[
            f"{TIME_COLUMN_NAME}_steadystate"
        ],
        (STEADY_STATE, data_column): steadystate[f"{data_column}_steadystate"],
        ("Fit", "params"): params,
        ("Fit", "R^2"): r_squared,
        ("Fit", "drift"): 0,
        ("Fit", 'good_fit'): r_squared > 0.95

    }

    if fit_function_name == "sigmoid_drift":
        # kinetics["Velocity", TIME_COLUMN_NAME] = pd.to_timedelta(params[2], unit="s")
        # kinetics["Velocity", data_column] = params[0]/2
        # kinetics["Velocity", "Max"] = params[1] *  params[0] /4 #max slope, but not k which is just steepness.
        alpha = 0.95
        kinetics[STEADY_STATE, data_column] = params[0] * alpha
        t_alpha = logistic_time_to_fraction(*params[0:3], alpha=alpha)
        kinetics[STEADY_STATE, TIME_COLUMN_NAME] = pd.to_timedelta(
            t_alpha, unit="h"
        )
        kinetics["Fit", "drift"] = params[3]
        # kinetics["Lag", TIME_COLUMN_NAME]= pd.to_timedelta(params[2], unit="s")
        # kinetics["Lag", data_column] = _sigmoid_drift(params[2],*params)

    return pd.Series(kinetics)


def kinetic_analysis(
    data: pd.DataFrame,
    group_by=["Name", "Read"],
    data_column=None,
    fit_function_name=DEFAULT_FIT_FUNCTION_NAME,
) -> pd.DataFrame:
    if data_column is None:
        data_column = get_active_data_column()

    if (
        data[READ_COLUMN_NAME].unique().size > 1
        and READ_COLUMN_NAME not in group_by
    ):
        log.warning(
            "Kinetic analysis is not grouped on `Read`, but multiple different read types exist in the data. This is probably not what you want."
        )

    group_by_with_wells = group_by.copy()
    if "Well" not in group_by:  # kinetics should be calculated per well, THEN averaged
        # print(f"WARNING: \"Well\" not explicit in group_by ({group_by}). Calculating individual replicate fits first, then averaging")
        group_by_with_wells += ["Well"]
        kinetics = kinetic_analysis(
            data,
            group_by=group_by_with_wells,
            fit_function_name=fit_function_name,
            data_column=data_column,
        )
        if "Well" not in group_by:
            print("PROVIDING AVERAGED KINETICS")
            return get_average_kinetics(kinetics.reset_index(), group_by=group_by, data_column=data_column)

    # Use apply with a lambda to pass group_keys explicitly
    kinetics = data.groupby(group_by, sort=False).apply(
        lambda g: kinetic_analysis_per_well(
            g,
            group_keys=g.name,  # Pass the group keys directly
            data_column=data_column,
            fit_function_name=fit_function_name,
        ),
        include_groups=False,
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
    group_by=["Name", READ_COLUMN_NAME],
    precision: float = 2,
):
    if kinetics is None:
        kinetics = kinetic_analysis(data, group_by=group_by)

    kinetics_styled = kinetics.style.format(precision=precision).format(
        _format_timedelta,
        subset=pd.IndexSlice[:, pd.IndexSlice[:, TIME_COLUMN_NAME]],
    )
    return kinetics_styled


def plot_kinetics_by_well(
    data: pd.DataFrame,
    kinetics: pd.DataFrame,
    group_by: list[str],
    x: str = TIME_COLUMN_NAME,
    y: str = None,
    show_data: bool = True,
    show_fit: bool = False,
    show_velocity: bool = False,
    show_mean: bool = False,
    annotate: bool = False,
    fit_function_name=DEFAULT_FIT_FUNCTION_NAME,
    **kwargs,
):
    """
    Typical usage:

    > tk = kinetic_analysis(data=data, data_column="BackgroundSubtracted")
    > g = sns.FacetGrid(tk, col="Well", col_wrap=2, sharey=False, height=4, aspect=1.5)
    > g.map_dataframe(plot_kinetics, show_fit=True, show_velocity=True)
    """
    if y is None:
        y = get_active_data_column()

    fit_function = FIT_FUNCTIONS[fit_function_name]
    log.debug(f"Plotting kinetics for group: {group_by}")

    colors = sns.color_palette("Set2")

    data_alpha = 0.25
    if show_data < 1:
        data_alpha = show_data  # Make data invisible. We still plot it to generate the axes.
    ax = sns.scatterplot(
        data=data, x=x, y=y, color=colors[2], alpha=data_alpha
    )

    # TODO: This silently just plots the first group_by
    data_index = data.iloc[0].loc[group_by]
    kinetics = kinetics.loc[*data_index]

    log.debug(f"Data index: {data_index.values}")

    if kinetics.isna().any():
        log.info(f"Kinetics information not available for {data_index}.")
        return

    if show_mean:
        log.debug("Plotting data mean")
        sns.scatterplot(
            data=_data_mean(data, data_col=y),
            x=x,
            y=y,
            color=colors[4],
            alpha=0.5,
            ax=ax,
        )

    # ax_ylim = (
    #     ax.get_ylim()
    # )  # Use this to run lines to bounds later, then restore them before returning.

    if show_fit:
        params = kinetics["Fit", "params"]
        sns.lineplot(
            x=data[x],
            y=fit_function(data[x].dt.total_seconds()/3600.0, *params),
            linestyle="--",
            color=colors[3],
            # alpha=0.5,
            ax=ax,
        )
    #     sns.lineplot(data=data, x=x, y=y, linestyle="--", c="red", alpha=0.5)

    # Max Velocity
    # maxV_x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    # The tangent line to the inflection point is:y= V_slope * (t-t0) + V_max_data
    maxV_y = (
        kinetics["Velocity", "Max"]
        * (
            data[x] - kinetics["Velocity", x]
        ).dt.total_seconds()/3600.0
        + kinetics["Velocity", y]
    )
    if fit_function_name == "sigmoid_drift":
        params = kinetics["Fit", "params"]
        maxV_y = maxV_y + params[3]*(data[x].dt.total_seconds()/3600.0 - params[4])

    sns.lineplot(
        x=data[x].loc[(maxV_y > 0) & (maxV_y < data[y].max())],
        y=maxV_y[(maxV_y > 0) & (maxV_y < data[y].max())],
        linestyle="--",
        color=colors[1],
        ax=ax,
    )

    maxV = kinetics["Velocity", "Max"]
    maxV_s = kinetics["Velocity", x]
    maxV_d = kinetics["Velocity", y]

    # Time to Steady State
    ss_s = kinetics[STEADY_STATE, x]
    if ss_s < max(data[x]): # Dont plot vert line if ss goes beyond reaction time
        ax.axvline(ss_s, c=colors[3], linestyle="--")
    ax.axhline(kinetics[STEADY_STATE, y], c=colors[3], linestyle="--")

    # # Range
    # ax.axhline(decile_upper, c=colors[7], linestyle="--")
    # ax.axhline(decile_lower, c=colors[7], linestyle="--")

    if annotate:
        # Plot the text annotations on the chart
        ax.annotate(
            f"$V_{{max}} =$ {maxV:.2f} u/h",
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
            timple.timedelta.timedelta2num(kinetics["Lag", x])
        )
        ax.annotate(
            f"$t_{{lag}} =$ {lag_label} h",
            (kinetics["Lag", x], kinetics["Lag", y]),
            xytext=(12, 0),
            textcoords="offset points",
            ha="left",
            va="center",
        )

        ss_label = f.format_data(
            timple.timedelta.timedelta2num(
                kinetics[STEADY_STATE, x]
            )
        )
        ax.annotate(
            f"$t_{{steady state}} =$ {ss_label} h",
            (
                min(kinetics[STEADY_STATE, x], max(data[x])),
                kinetics[STEADY_STATE, y],
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
        velocity_ax.set_ylabel("$V (u/h)$")
        velocity_ax.set_ylim((0, velocity[y].max() * 2))

    # ax.set_ylim(ax_ylim)

    _plot_timedelta(ax)

def get_average_kinetics(kinetics, group_by, data_column=None, time=TIME_COLUMN_NAME):
    if data_column is None:
        data_column = get_active_data_column()

    # These are the columns to average
    cols_to_avg = [
        ('Velocity', time),
        ('Velocity', data_column),
        ('Velocity', 'Max'),
        ('Lag', time),
        ('Lag', data_column),
        ('Steady State', time),
        ('Steady State', data_column),
        ('Fit', 'params'),
        ('Fit', 'R^2'),
        ('Fit', 'drift')
    ]
    avg_df = kinetics.groupby(group_by)[cols_to_avg].mean() #.reset_index()
    return avg_df

def plot_kinetics(
    data: pd.DataFrame,
    group_by=["Name", READ_COLUMN_NAME],
    kinetics: pd.DataFrame = None,
    show_data: bool = True,
    show_fit: bool = True,
    show_velocity: bool = False,
    show_mean: bool = True,
    annotate: bool = True,
    fit_function_name: str = DEFAULT_FIT_FUNCTION_NAME,
    data_column: str = None,
    **kwargs,
):
    if data_column is None:
        data_column = get_active_data_column()

    if kinetics is None: #TODO: force that the group_bys are the same if input is provided.
        kinetics = kinetic_analysis(
            data,
            group_by=group_by,
            fit_function_name=fit_function_name,
            data_column=data_column,
        )

    if "col_wrap" not in kwargs:
        kwargs["col_wrap"] = 3

    # TODO: this just plots first gain without issuing a warning.
    g = sns.FacetGrid(data, col=group_by[0], height=4, aspect=1.5, **kwargs)
    g.map_dataframe(
        plot_kinetics_by_well,
        y=data_column,
        kinetics=kinetics,
        group_by=group_by,
        show_data=show_data,
        show_fit=show_fit,
        show_velocity=show_velocity,
        show_mean=show_mean,
        annotate=annotate,
        fit_function_name=fit_function_name,
    )
    if "_norm" in data_column:
        g.set_ylabels("Normalized Fluorescence (RFU)")
    else:
        g.set_ylabels("Fluorescence (AFU)")
    return g, kinetics


def plot_steadystate(
    data: pd.DataFrame,
    x="Name",
    show_points=True,
    show_plot=False,
    fit_function_name=DEFAULT_FIT_FUNCTION_NAME,
    data_column=None,
    **kwargs,
):
    if data_column is None:
        data_column = get_active_data_column()

    steady_state = find_steady_state(
        data,
        group_by=[
            "Well",
            READ_COLUMN_NAME,
            x,
        ],  # if we put in x here, no need to merge as below.m
        fit_function_name=fit_function_name,
        data_column = data_column
    ).reset_index()
    ss_name = f"{data_column}_steadystate"

    # Merge in the original data to get access to all the platemap columns,
    # but drop duplicates so we don't have the full table with every single data point, but
    # rather one row per Well (which will have the steady state for that well attached).
    # data_with_steady_state = steady_state.merge(
    #     data, on=["Well", GAIN_COLUMN_NAME], how="left"
    # )
    # data_with_steady_state = data_with_steady_state.drop_duplicates(
    #     subset=["Well", GAIN_COLUMN_NAME]
    # )

    if "sharex" not in kwargs:
        kwargs["sharex"] = False

    if "col_wrap" not in kwargs and "col" in kwargs:
        kwargs["col_wrap"] = 2

    g = sns.catplot(
        data=steady_state,
        x=x,
        y=ss_name,
        kind="bar",
        height=4,
        aspect=1.5,
        **kwargs,
    )

    if show_points:
        # Pull hue directly, because we don't want to drop in the full set of facet args that
        # might have been provided.
        g.map_dataframe(
            sns.stripplot,
            x=x,
            y=ss_name,
            hue=kwargs["hue"] if "hue" in kwargs else None,
            palette=sns.color_palette(
                desat=0.5
            ),  # Set palette directly because this plot gets transformed hue data which isn't categorical.
            size=8,
            linewidth=2,
            marker=".",
            jitter=False,
            dodge=False,
        )

    g.set_xticklabels(rotation=90)
    g.set_ylabels("Steady State Fluorescence (RFU)")

    if show_plot:
        plt.tight_layout()
        plt.show()

    return g

def get_avg_ctrl(ctrls, min_to_avg = 60, time_int = 5):
    """
    Compute the average control intensity over the final portion of the time trace.

    Parameters
    ----------
    ctrls : pandas.DataFrame
        Control measurements for a single (Name, Read) group.
        Must contain a column named "Data" and be time-ordered.
    min_to_avg : int, optional
        Number of minutes from the end of the trace to include in the average.
        Default is 60 minutes.
    time_int : int, optional
        Time interval (in minutes) between measurements.
        Default is 5 minutes.

    Returns
    -------
    float
        The average value of the final `min_to_avg` minutes of control data.

    Notes
    -----
    This function assumes rows are evenly spaced in time.
    """
    num_pts = round(min_to_avg / time_int)
    return np.mean(ctrls[::-num_pts]["Data"])

def split_data_ctrls(df, ctrl_names = 'HPTS (10 uM)'):
    """
    Split a DataFrame into control rows and non-control rows.

    Parameters
    ----------
    df : pandas.DataFrame
        Full experimental dataset, containing a column named "Name".
    ctrl_names : str or list of str, optional
        Name(s) identifying control samples. Default is 'HPTS (10 uM)'.

    Returns
    -------
    data : pandas.DataFrame
        All rows not matching the control names.
    ctrls : pandas.DataFrame
        All rows matching the control names.
    """
    ctrls = df[df.Name.isin([ctrl_names])]
    data = df[~df.Name.isin([ctrl_names])]
    return data, ctrls

def normalize_data_to_controls(data:pd.DataFrame,
                               ctrl_name: str = 'HPTS (10 uM)',
                               output_col: str = "data_normalized",
                               set_active: bool = True) -> pd.DataFrame:
    """
    Normalize experimental data to control intensities.

    The function:
        1. Splits the data into control rows and non-control rows.
        2. Computes an average control intensity for each (Name, Read) group.
        3. Merges control averages back into the experimental data.
        4. Creates a new column: ``data_norm`` = ``Data`` / control average.

    Parameters
    ----------
    data : pandas.DataFrame
        Full dataset with a "Data" column and at least "Name" and "Read".
    ctrl_name : str, optional
        Name of the control sample used for normalization.

    Returns
    -------
    pandas.DataFrame
        A copy of the input data with an additional column:
        - ``data_norm`` : normalized intensity values.

    Notes
    -----
    After normalization, pass ``data_column='data_norm'`` to downstream plotting/
    analysis functions to use normalized data.
    """
    data, ctrls = split_data_ctrls(data, ctrl_names=ctrl_name)

    ctrl_avgs = ctrls.groupby(["Name", "Read"]).apply(get_avg_ctrl, include_groups=False).rename(
        "ctrl_ints").reset_index()
    df_merge = data.merge(ctrl_avgs[["Read", "ctrl_ints"]], on=["Read"], how="left")
    df_merge[output_col] = df_merge["Data"] / df_merge["ctrl_ints"]

    if set_active:
        set_active_data_column(output_col)

    print(f"Data Normalized to {ctrl_name} in col {output_col}. "
          f"The active column for subsequent operations is: {get_active_data_column()}")

    return df_merge

##### PLOTTING METHODS ####

def plot_summary(
    data: pd.DataFrame,
    experiment_split="Name",
    show_plot=True,
    data_column=None,
    fit_function_name=DEFAULT_FIT_FUNCTION_NAME,
):
    """
    Generate a summary visualization of plate reader data across experiments.

    This function produces:
        1. Line plots of fluorescence (or specified data column) over time for each experiment.
        2. Bar plots summarizing kinetic parameters (steady-state, max velocity, and optionally drift)
           for each experiment and group, split by the specified column (default "Name").

    Args:
        data (pd.DataFrame): Input plate reader data.
        experiment_split (str): Column name used to split lines and bars by group (default "Name").
        show_plot (bool): Whether to show plot explicitly
        data_col (str): Column in `data` to plot as y-axis (default ACTIVE_DATA_COLUMN_NAME).
        fit_function_name (str): Fit type for kinetics ("sigmoid" or "sigmoid_drift").

    Returns:
        None: Plots are displayed and optionally arranged in subplots.
    """
    if data_column is None:
        data_column = get_active_data_column()

    # --- Identify experiments ---
    experiments = [("Experiment", data)]
    if "Experiment" in data:
        experiments = list()
        for experiment in data["Experiment"].unique():
            experiments.append(
                (experiment, data[data["Experiment"] == experiment])
            )

    exp_count = len(experiments)

    # --- Configure number of columns based on fit type ---
    ncols = {"sigmoid_drift": 4, "sigmoid": 3}
    fig, axes = plt.subplots(
        exp_count,
        ncols[fit_function_name],
        figsize=(15, 5 * exp_count),
        squeeze=False,
    )

    # --- Loop over each experiment and create plots ---
    for i, (experiment, df) in enumerate(experiments):
        ax_c = sns.lineplot(
            data=df,
            x=TIME_COLUMN_NAME,
            y=data_column,
            hue=experiment_split,
            legend=False,
            ax=axes[i, 0],
        )
        _plot_timedelta(axes[i, 0])
        ax_c.set_title(experiment)
        ax_c.set_xlabel(TIME_COLUMN_NAME)
        if "_norm" in data_column:
            ax_c.set_ylabels("Normalized Fluorescence (RFU)")
        else:
            ax_c.set_ylabels("Fluorescence (AFU)")


        # Prepare kinetic summary data
        ys_to_plot = [(STEADY_STATE, data_column), ("Velocity", "Max")]
        if fit_function_name == "sigmoid_drift":
            ys_to_plot.append(("Fit", "drift"))

        ss = kinetic_analysis(
            df,
            group_by=[experiment_split, "Well", READ_COLUMN_NAME],
            fit_function_name=fit_function_name,
            data_column=data_column,
        )

        # Helper function to plot bar plots with strip plots overlaid
        def plot_fit_params(ss, x, y, hue, ax, label=""):
            ax_ss = sns.barplot(
                data=ss,
                x=x,
                y=y,
                hue=hue,
                ax=ax,
            )
            palette = {True: "#1f77b4", False: "#ff7f0e"}
            sns.stripplot(x=x, y=y, data=ss, ax=ax, hue=("Fit", "good_fit"), palette=palette)
            ax_ss.tick_params(axis="x", labelsize="x-small", rotation=90)
            ax_ss.set_xlabel(x)
            ax_ss.set_ylabel(y)
            ax_ss.set_title(f"{y}, {label}")

        # Plot kinetics parameters in remaining columns
        for ax, cur_y in zip(axes[i, 1:], ys_to_plot):
            plot_fit_params(
                ss=ss,
                x=experiment_split,
                y=cur_y,
                hue=experiment_split,
                ax=ax,
                label=fit_function_name,
            )

        if show_plot:
            plt.tight_layout()
            plt.show()

def plot_main_effects(data, factors, response_var):
    """
    Plot main effects of experimental factors on a response variable.
    Parameters
    ----------
    data : pandas.DataFrame
    factors : list of str
        Names of columns in `data` to be treated as experimental factors.
        A separate subplot is created for each factor and the response is
        aggregated over all other dimensions.
    response_var : str
        `data` containing the response variable
        (e.g. steady-state signal, Vmax, yield).
    """
    n = len(factors)

    # Determine subplot grid size
    cols = 3
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(12, 4 * rows))
    axes = axes.flatten()  # flatten in case of 1 row

    for ax, factor in zip(axes, factors):
        sns.pointplot(
            x=factor,
            y=response_var,
            data=data,
            errorbar="sd",
            markers="o",
            linestyles="-",
            ax=ax
        )
        ax.set_title(f"Main Effect of {factor} on {response_var}")
        ax.set_xlabel(factor)
        ax.set_ylabel(response_var)

    # Hide any empty subplots
    for i in range(len(factors), len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    plt.show()

def plot_all_interactions(data, factors, response_var, show_all_points=True):
    """
    Plot pairwise interaction effects between experimental factors.

    Parameters
    ----------
    data : pandas.DataFrame
        Must include all columns listed in `factors` and the `response_var` column.
    factors : list of str
        Names of columns in `data` to be treated as experimental factors.
        All pairwise combinations of these factors are plotted as interactions.
    response_var : str
        `data` column containing the response variable
        (e.g. steady-state signal, Vmax, yield).
    show_all_points : bool, default True
        If True, plot all individual data points as scatter markers in addition
        to the mean ± SD lines. If False, only the smoothed/aggregated lines
        with error bars are shown, which can be clearer for dense datasets.
    """
    factor_pairs = list(combinations(factors, 2))
    n = len(factor_pairs)

    cols = 3
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(12, 5 * rows))
    axes = axes.flatten()

    for ax, (factor_x, factor_hue) in zip(axes, factor_pairs):

        levels = np.sort(data[factor_hue].unique())

        # Build a continuous colormap (viridis is great for continuous data)
        cmap = sns.color_palette("viridis", as_cmap=True)

        # Normalize hue values to [0,1] for colormap mapping
        if np.issubdtype(data[factor_hue].dtype, np.number):
            norm_levels = (levels - levels.min()) / (levels.max() - levels.min() + 1e-9)
        else:
            # For categorical: evenly spaced colors
            norm_levels = np.linspace(0, 1, len(levels))

        color_map = {level: cmap(norm) for level, norm in zip(levels, norm_levels)}

        # Plot per hue level
        for level in levels:
            subset = data[data[factor_hue] == level]

            if show_all_points:
                # Scatter points
                sns.scatterplot(
                    data=subset,
                    x=factor_x,
                    y=response_var,
                    ax=ax,
                    color=color_map[level],
                    label=f"{factor_hue}: {level}"
                )

            # CI-smoothed line, same color
            sns.lineplot(
                data=subset,
                x=factor_x,
                y=response_var,
                ax=ax,
                color=color_map[level],
                errorbar="sd"
            )

        ax.set_title(f"Interaction: {factor_x} × {factor_hue}")
        ax.set_xlabel(factor_x)
        ax.set_ylabel(response_var)

    # Hide unused axes
    for i in range(len(factor_pairs), len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    plt.show()


def plot_param_heatmaps(
    df: pd.DataFrame,
    row_col: str,
    col_col: str,
    facet_col: str,
    value_col,
    title: str | None = None,
    cmap: str = "viridis",
    annot: bool = True,
) -> None:
    """
    Plot a series of heatmaps for an n-dimensional parameter sweep with averaged replicates.

    This function expects a dataframe where each row corresponds to a single
    measurement from a combinatorial experiment (e.g. Mg, Ribo, T7RNAP sweeps).
    It will:
      1. Group by (row_col, col_col, facet_col) and average replicate measurements
         in value_col.
      2. Compute a global minimum and maximum across all facets so that every
         heatmap shares the same color scale.
      3. Create one heatmap per unique value of facet_col, arranged in a single row.

    Parameters
    ----------
    df : pandas.DataFrame
        Long-form dataframe containing the experimental data. Must include
        the columns specified in row_col, col_col, facet_col, and value_col.
    row_col : str
        Column name to use as the heatmap y-axis (heatmap rows).
        Example: "Ribo_uM".
    col_col : str
        Column name to use as the heatmap x-axis (heatmap columns).
        Example: "Mg_mM".
    facet_col : str
        Column name defining the different “slices” of the experiment.
        One heatmap will be drawn for each unique value in this column.
        Example: "T7poly_ng/ul".
    value_col : str
        Column name containing the measurement to visualize in the heatmap
        (e.g. steady-state value, Vmax, area under curve).
    title : str, optional
        Figure-level title to display above all heatmaps. If None, no
        super-title is added.
    cmap : str, default "viridis"
        Name of the matplotlib colormap to use for the heatmaps.
    annot : bool, default True
        If True, write the numeric value into each cell of the heatmap.

    Returns
    -------
    None
        The function creates a matplotlib figure and displays it with plt.show().
        It does not return any objects.

    Notes
    -----
    - Replicates are averaged using a groupby over (row_col, col_col, facet_col).
    - All facets share a single colorbar and identical color scaling, which makes
      visual comparison across slices straightforward.
    """
    if type(value_col) is tuple:
        df["data_column"] = df[value_col] # resets just in case is multiindex
        value_col = "data_column"
    # Average replicates over the parameter combinations
    means = df.groupby([row_col, col_col, facet_col], as_index=False)[value_col].mean()

    # Shared color scale across all facets
    data_min = means[value_col].min()
    data_max = means[value_col].max()

    facet_values = sorted(means[facet_col].unique())
    n_plots = len(facet_values)

    fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 5), squeeze=False)

    # Single shared colorbar axis
    cbar_ax = plt.gcf().add_axes([0.93, 0.25, 0.02, 0.5]) if n_plots > 0 else None

    for idx, fval in enumerate(facet_values):
        slice_df = means[means[facet_col] == fval]

        # Pivot to 2D matrix for heatmap
        pivot = slice_df.pivot(index=row_col, columns=col_col, values=value_col)

        ax = axes[0, idx]
        sns.heatmap(
            pivot,
            annot=annot,
            cmap=cmap,
            vmin=data_min,
            vmax=data_max,
            cbar=(idx == n_plots - 1),
            cbar_ax=cbar_ax if idx == n_plots - 1 else None,
            ax=ax,
        )
        ax.set_title(f"{facet_col} = {fval}")
        ax.set_xlabel(col_col)
        ax.set_ylabel(row_col)

    if title:
        plt.suptitle(title)

    plt.tight_layout(rect=[0, 0, 0.9, 0.95])
    plt.show()

# import statsmodels.api as sm # not in current dependencies
# from statsmodels.formula.api import ols # not in current dependencies
# def run_full_factorial_anova(
#     df,
#     factors,
#     response_var,
#     max_interaction_order=None,
#     rename_conflicts=True,
# ):
#     """
#     Fit an OLS model with all main effects and interaction terms up to a given order,
#     then run an ANOVA table.
#
#     Parameters
#     ----------
#     df : pandas.DataFrame
#         Input data frame containing the response and factor columns.
#     factors : list of str
#         Column names to treat as factors (can contain characters like '/' etc.).
#     response_var : str
#         Column name of the response variable (e.g. "steady_state").
#     max_interaction_order : int or None, default None
#         Maximum interaction order to include.
#         - If None, use len(factors) (full factorial, including k-way interaction).
#         - If 2, include all main effects and all 2-way interactions only.
#     rename_conflicts : bool, default True
#         If True, create “safe” versions of factor names for the formula
#         (replacing non-alphanumeric chars with '_') and work on a copy of df.
#
#     Returns
#     -------
#     model : statsmodels.regression.linear_model.RegressionResultsWrapper
#         Fitted OLS model.
#     anova : pandas.DataFrame
#         ANOVA table from statsmodels.stats.anova_lm().
#     """
#     if type(response_var) is tuple:
#         df["data_column"] = df[response_var] # resets just in case is multiindex
#         response_var = "data_column"
#
#     if max_interaction_order is None:
#         max_interaction_order = len(factors)
#
#     # Work on a copy to avoid mutating the original df
#     data = df.copy()
#
#     # Build mapping from original factor names to formula-safe names
#     safe_name_map = {}
#     for f in factors:
#         if rename_conflicts:
#             safe = re.sub(r'[^0-9a-zA-Z_]', '_', f)
#         else:
#             safe = f
#         # Avoid collisions if two columns sanitize to same name
#         if safe in safe_name_map.values():
#             # Append a numeric suffix until unique
#             base = safe
#             i = 1
#             while f"{base}_{i}" in safe_name_map.values():
#                 i += 1
#             safe = f"{base}_{i}"
#         safe_name_map[f] = safe
#
#     # Apply renaming in the working copy
#     data = data.rename(columns=safe_name_map)
#
#     # Response name is assumed to already be formula-safe;
#     # if not, you can apply the same cleaning logic.
#     response_safe = re.sub(r'[^0-9a-zA-Z_]', '_', response_var)
#
#     # Build formula parts
#     safe_factors = [safe_name_map[f] for f in factors]
#
#     # Main effects
#     terms = safe_factors.copy()
#
#     # Interaction terms up to max_interaction_order
#     for order in range(2, max_interaction_order + 1):
#         for combo in combinations(safe_factors, order):
#             terms.append(":".join(combo))
#
#     rhs = " + ".join(terms)
#     formula = f"{response_safe} ~ {rhs}"
#
#     # Fit model and run ANOVA
#     model = ols(formula, data=data).fit()
#     anova = sm.stats.anova_lm(model)
#
#     return model, anova

def compute_standard_curve(
    data: pd.DataFrame, sc_type="STD", include_mean=True
):
    std = data[data["Type"] == sc_type]

    def curve(group):
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            group["Conc/Dil"], group[get_active_data_column()]
        )
        return pd.Series(
            {
                "slope": slope,
                "intercept": intercept,
                "min": group["Conc/Dil"].min(),
                "max": group["Conc/Dil"].max(),
                "R": r_value,
                "R^2": r_value**2,
            }
        )

    curves = (
        std.groupby(["Type", READ_COLUMN_NAME])
        .apply(curve, include_groups=False)
        .reset_index()
    )

    if include_mean:
        mean = (
            std.groupby(["Type"])
            .apply(curve, include_groups=False)
            .reset_index()
        )
        mean[READ_COLUMN_NAME] = "Mean"
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
        on=[READ_COLUMN_NAME],
        how="left",
        suffixes=("", "_y"),
    )
    conc["Concentration"] = (
        conc[get_active_data_column()] - conc["intercept"]
    ) / conc["slope"]

    # Exclude data points which fall outside the standard curve.
    conc["In SC"] = False
    conc.loc[
        (conc["Concentration"] >= conc["min"])
        & (conc["Concentration"] <= conc["max"]),
        "In SC",
    ] = True

    if dilution_col in data.columns:
        conc["Original Concentration"] = (
            conc["Concentration"] * conc[dilution_col]
        )

    return conc


def plot_standard_curve(data: pd.DataFrame, sc_type="STD", **kwargs):
    curves = compute_standard_curve(data, sc_type)

    std_data = data.merge(curves, on=["Type", READ_COLUMN_NAME])
    std_data["StdCurvePoint"] = (
        std_data["Conc/Dil"] * std_data["slope"] + std_data["intercept"]
    )

    if "col" not in kwargs:
        kwargs["col"] = READ_COLUMN_NAME

    if "hue" not in kwargs:
        kwargs["hue"] = READ_COLUMN_NAME

    g = sns.FacetGrid(data=std_data, **kwargs)
    g.map(sns.scatterplot, "Conc/Dil", get_active_data_column())
    g.map(sns.lineplot, "Conc/Dil", "StdCurvePoint", linestyle="--")
    g.set_titles(col_template="{col_name}")
    g.add_legend()

    def annotate(data, **kwargs):
        r2 = data["R^2"].iloc[0]
        plt.text(0.7, 0.1, f"$R^2 = {r2:.3f}$", transform=plt.gca().transAxes)

    g.map_dataframe(annotate)

    unit = ""
    if "Unit" in std_data:
        unit = std_data["Unit"].iloc[0]
        unit = f"({unit})"

    g.set_xlabels(f"Concentration {unit}")
    g.set_ylabels("Relative Fluorescence Units (RFU)")

    return g


def plot_concentration(
    data: pd.DataFrame,
    x="Name",
    y=None,
    hue=READ_COLUMN_NAME,
    show_outliers=False,
    show_box=True,
    log_scale=True,
    **kwargs,
):
    conc = compute_concentration(data)

    y_column = y if y is not None else "Concentration"
    if "Original Concentration" in conc.columns:
        y_column = "Original Concentration"

    def plot_with_in_sc(
        data, x, y, color, label, show_outliers, show_box, log_scale
    ):
        ax = sns.stripplot(
            data=data[data["In SC"]],
            x=x,
            y=y,
            color=color,
            label=label,
            log_scale=log_scale,
            dodge=True,
            jitter=False,
        )

        # means = data[data["In SC"]].groupby(["Name"])["Concentration"].mean().reset_index()
        if show_box:
            sns.boxplot(
                data=data[data["In SC"]],
                x=x,
                y=y,
                # meanline=True,
                # showmeans=True,
                showcaps=False,
                showbox=True,
                showfliers=True,
                meanprops={"ls": "-", "lw": 2},
                medianprops={"visible": True},
                whiskerprops={"alpha": 0.25},
                boxprops={"alpha": 0.5},
                flierprops={"alpha": 0.5},
                # fliersize=50,
                zorder=10,
                fill=False,
                legend=True,
                dodge=True,
                width=0.4,
                color=color,
                ax=ax,
            )

        if show_outliers:
            sns.stripplot(
                data=data[~data["In SC"]],
                x=x,
                y=y,
                edgecolor="auto",
                linewidth=1,
                alpha=0.4,
                color=color,
                ax=ax,
                dodge=True,
                jitter=False,
            )

    g = sns.FacetGrid(data=conc, hue=hue, **kwargs)
    g.map_dataframe(
        plot_with_in_sc,
        x=x,
        y=y_column,
        show_outliers=show_outliers,
        show_box=show_box,
        log_scale=log_scale,
    )
    g.tight_layout()
    g.add_legend()

    # g = sns.catplot(data=conc[conc["In SC"]], x=x, y=y_column, hue=hue)
    # g.map_dataframe(sns.stripplot, data=conc[~conc["In SC"]], x=x, y=y_column, color="red")

    # ax = sns.stripplot(data=conc[conc["In SC"]], x=x, y=y_column, hue=hue, **kwargs)
    # sns.stripplot(data=conc[~conc["In SC"]], x=x, y=y_column, hue=hue, **kwargs)

    # means = conc.groupby(["Name"])["Concentration"].mean().reset_index()

    # sns.boxplot(
    #     data=means,
    #     x="Name",
    #     y="Concentration",
    #     showmeans=True,
    #     meanline=True,
    #     meanprops={"color": "k", "ls": "--", "lw": 1, "alpha": 0.5},
    #     medianprops={"visible": False},
    #     whiskerprops={"visible": False},
    #     zorder=10,
    #     showfliers=False,
    #     showbox=False,
    #     showcaps=False,
    #     fill=False,
    #     legend=True,
    #     dodge=False,
    #     # ax=ax
    # )

    # if "Unit" in conc:
    #     ax.set_ylabel(f"Concentration ({conc['Unit'].iloc[0]})")

    return g


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
            df[col] = df[col].dt.total_seconds()/3600.0
            df = df.rename(columns={col: f"{col} (h)"})

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
        p[TIME_COLUMN_NAME] = p["Clock Time"] - start_time
        plate_data.append(p)

    merged_data = pd.concat(plate_data)
    merged_data["Plate"] = plates[0]

    return merged_data
