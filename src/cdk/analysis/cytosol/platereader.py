"""
Plate Reader Module

This module provides support for loading and analyzing data from various plate readers. Currently supported are:
+ BioTek Cytation 5
+ Revvity Envision Nexus
+ Promega Glomax Discover

Usage:
*tbd*
"""

import io
import re
import os.path
from pathlib import Path
import logging
from typing import Union, Optional, NamedTuple
from enum import Enum, auto

import numpy as np
import pandas as pd

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
    LOWEST_VELOCITY = auto(),
    MAXIMUM_VALUE = auto(),
    VELOCITY_INTERCEPT = auto()


class PlateReaderData(NamedTuple):
    """
    A named tuple representing plate reader data and associated plate map information.
    """

    data: pd.DataFrame
    platemap: Optional[pd.DataFrame] = None


_timple = timple.Timple()


def load_platereader_data(
    data_file: DataFile, platemap_file: Optional[DataFile] = None
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
    filename = os.path.basename(data_file).lower()

    if filename.startswith("cytation"):
        data = read_cytation(data_file)
    elif filename.startswith("envision"):
        data = read_envision(data_file)
    # elif filename_lower.startswith("glomax"):
    #     return read_glomax(os.path.dirname(data_file))
    else:
        raise ValueError(f"Unsupported plate reader data file: {data_file}")

    platemap = None
    if platemap_file is not None:
        platemap = read_platemap(platemap_file)
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
    platemap = platemap.convert_dtypes(convert_floating=False)
    
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


def read_cytation(data_file: DataFile, sep="\t") -> pd.DataFrame:
    log.debug(f"Reading Cytation data from {data_file}")
    # read data file as long string
    data = ""
    with open(data_file, "r", encoding="latin1") as file:
        data = file.read()

    # extract indices for Proc Details, Layout
    procidx = re.search(r"Procedure Details", data)
    layoutidx = re.search(r"Layout", data)
    readidx = re.search(r"^(Read\s)?\d+,\d+", data, re.MULTILINE)

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

    sep = r"(?:Read\s\d+:)?(?:\s\d{3},\d{3}(?:\[\d\])?)?" + sep

    for readidx in re.finditer(r"^(Read\s)?\d+,\d+.*\n", data, re.MULTILINE):
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
        r["Temperature (C)"] = r["T°"].str.extract(r"(\d+)").astype(float)
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
        lambda well: f"{well['Row']}:{well['Column']}", axis=1
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


def plot_setup() -> None:
    _timple.enable()
    pd.set_option('display.float_format', '{:.2f}'.format)


def _plot_timedelta(plot: sns.FacetGrid | mpl.axes.Axes) -> sns.FacetGrid:
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
        data=data, x="Time", y="Data", row="Row", col="Column", kind="line"
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
    if "Experiment" in data.columns and by_experiment:
        kwargs["row"] = "Experiment"

    g = plot_curves(
        data=data, x="Time", y="Data", hue="Name", col="Read", **kwargs
    )

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
    g = sns.relplot(data=data, x=x, y=y, hue=hue, kind="line", **kwargs)
    _plot_timedelta(g)

    x_label, y_label = labels
    if x_label:
        g.set_xlabels(x_label)
    if y_label:
        g.set_ylabels(y_label)

    # Set simple row and column titles, if we're faceting on row or column.
    # The join means the punctuation only gets added if we have both.
    row_title = "{row_name}" if "row" in kwargs else ""
    col_title = "{col_name}" if "col" in kwargs else ""
    g.set_titles(": ".join(filter(None, [row_title, col_title])))

    return g


###
# Kinetics Analysis
# TODO: Perhaps split this out into a submodule.
###


def find_steady_state_for_well(well):
    well = well.sort_values("Time")
    pct_change = well["Data"].rolling(window=3).mean().pct_change()
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


def kinetic_analysis_per_well(
    data: pd.DataFrame, data_column="Data"
) -> pd.DataFrame:

    steadystate = find_steady_state_for_well(data)

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
    log.debug(f"{data['Well'].iloc[0]} Fitted params: {params}")

    # calculate velocities and velocity params
    v = data[data_column].diff() / data["Time"].dt.total_seconds().diff()

    maxV = v.max()
    maxV_d = data.loc[v.idxmax(), data_column]
    maxV_time = data.loc[v.idxmax(), "Time"]

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
    }

    return pd.Series(kinetics)
    # return kinetics


def kinetic_analysis(data: pd.DataFrame, data_column="Data") -> pd.DataFrame:
    kinetics = data.groupby(["Well", "Read"]).apply(
        functools.partial(kinetic_analysis_per_well, data_column=data_column)
    )
    return kinetics


def kinetic_analysis_summary(
    data: pd.DataFrame,
    data_column="Data",
    time_cutoff: int = 12000,
    label_order: list[str] = None,
    norm_label: str = None,
):
    def per_well_cleanup(df):
        cols = df.columns
        return df[["Well"] + list(cols[27:])].aggregate(lambda x: x.iloc[0])

    tk = kinetic_analysis(
        data=data, data_column=data_column, time_cutoff=time_cutoff
    )
    out = tk.groupby("Well").apply(per_well_cleanup).reset_index(drop=True)

    if label_order:
        out = out.set_index("Well").reindex(label_order).reset_index()

    # normalize max value (calculated by kinetics) if norm_label given
    if norm_label:
        norm = out[out["Well"] == norm_label][f"{data_column}_high_d"].values
        out["Normalized [%]"] = out[f"{data_column}_high_d"] / norm

    return out


def plot_kinetics_by_well(
    data: pd.DataFrame,
    kinetics: pd.DataFrame,
    x: str = "Time",
    y: str = "Data",
    show_fit: bool = False,
    show_velocity: bool = False,
    annotate: bool = False,
    **kwargs,
):
    """
    Typical usage:

    > tk = kinetic_analysis(data=data, data_column="BackgroundSubtracted")
    > g = sns.FacetGrid(tk, col="Well", col_wrap=2, sharey=False, height=4, aspect=1.5)
    > g.map_dataframe(plot_kinetics, show_fit=True, show_velocity=True)
    """
    colors = sns.color_palette("Set2")

    ax = sns.scatterplot(data=data, x=x, y=y, color=colors[2], alpha=0.5)

    well = data["Well"].iloc[0]
    read = data["Read"].iloc[0]
    kinetics = kinetics.loc[well, read]
    if (kinetics.isna()).any():
        log.info(f"Kinetics information not available for {well}.")
        return

    ax_ylim = (
        ax.get_ylim()
    )  # Use this to run lines to bounds later, then restore them before returning.

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

    ax.set_ylim(ax_ylim)

    _plot_timedelta(ax)


def plot_kinetics(data: pd.DataFrame, kinetics: pd.DataFrame, **kwargs):
    g = sns.FacetGrid(
        data, col="Well", col_wrap=2, sharey=True, height=4, aspect=1.5
    )
    g.map_dataframe(
        plot_kinetics_by_well,
        kinetics=kinetics,
        show_fit=True,
        show_velocity=False,
        annotate=True,
    )


def plot_steadystate(data: pd.DataFrame, hue="Experiment", **kwargs):
    steady_state = find_steady_state(data).reset_index()
    data_with_steady_state = data.merge(steady_state, on="Well", how="left")
    g = sns.catplot(
        data=data_with_steady_state,
        x="Name",
        y="Data_steadystate",
        hue=hue,
        kind="bar",
        col="Experiment",
        col_wrap=2,
        height=4,
        aspect=1.5,
        sharex=False,
        **kwargs
    )

    g.set_xticklabels(rotation=90)
    g.set_ylabels("Steady State Fluorescence (RFU)")
    return g


def export():
    # TODO: Write me
    pass