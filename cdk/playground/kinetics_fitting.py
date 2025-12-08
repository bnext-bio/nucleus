import pandas as pd
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import seaborn as sns
import numpy as np
import os
from cdk.analysis.cytosol import platereader as pr

# Testing out how we might fit to PURE timeseries data
# Load in data from experiments

def clean_gain_naming(df):
    df["gain"] = df.gain.apply(lambda x: x.replace(r"-485,528", ""))
    df["gain"] = df.gain.apply(lambda x: x.replace(r"GFP-F-", "GFP-"))
    df["gain"] = df.gain.apply(lambda x: x.replace(r"GFP-M-", "GFP-"))
    return df


def load_files():
    root_folder = "/Users/sharonnewman/Documents/general/data/platereader/calibration/"
    data_cyt5_name = os.path.join(root_folder, "20250908-112119-cytation5-pure-timecourse-gfp-exp5.txt")
    platemap_file = os.path.join(root_folder, "ot-pure-platemap-exp2.csv")

    data, platemap = pr.load_platereader_data(data_cyt5_name, platemap_file, platereader="biotek-cdk")
    data["gain"] = data["Read"]
    data = clean_gain_naming(data)
    ctrls = data[data.Name.isin(["HPTS_ctrls"])]
    data = data[~data.Name.isin(["HPTS_ctrls"])]
    return data, ctrls


def plot_kinetics(
        data: pd.DataFrame,
        params: list[float],
        x: str = "Time",
        y: str = "Data",
        show_data: bool = True,
        fit_function=pr._sigmoid,
        show_fit: bool = False,
        title: str = "",
):
    colors = sns.color_palette("Set2")

    data_alpha = 0.25
    if not show_data:
        data_alpha = show_data  # Make data invisible. We still plot it to generate the axes.
    ax = sns.scatterplot(
        data=data, x=x, y=y, color=colors[2], alpha=data_alpha
    )

    time_seconds = data["Time"].dt.total_seconds().to_list()
    est_fluorescence = fit_function(time_seconds, *params)

    sns.lineplot(
        x=data["Time"].to_list(),
        y=est_fluorescence,
        linestyle="--",
        color=colors[3],
        # alpha=0.5,
        ax=ax,
    )
    ax.axhline(params[0], c=colors[3], linestyle="--")
    plt.title(title)
    plt.show()


def _gompertz(t, L, k, t0):
    return L * np.exp(-np.exp(k * (t0 - t)))


def _gompertz_drift(t, L, k, t0, b, tau):
    """
    :param t:
    :param L:
    :param k:
    :param t0:
    :param b: drift rate
    :param tau: drift at onset
    :return:
    """
    return _gompertz(t, L, k, t0) + b * np.maximum(0, (t - tau))


def _weibull(t, A, tau, beta):
    return A * (1 - np.exp(-(t / tau) ** beta))


def _weibull_drift(t, L, k, t0, b, tau):
    """
    :param t:
    :param L:
    :param k:
    :param t0:
    :param b: drift rate
    :param tau: drift at onset
    :return:
    """
    return _weibull(t, L, k, t0) + b * np.maximum(0, (t - tau))


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
    return pr._sigmoid(t, L, k, t0) + b * (t - tau)


def fit_to_logistic(
        data: pd.DataFrame, data_column="Data", function=pr._sigmoid
):
    time = data["Time"].dt.total_seconds()
    signal_max_initial = np.max(data[data_column])
    k_initial = 0.007  #data_mean_v.max() / timestep / L_initial * 4
    t0_initial = 6000  #data["Time"].dt.total_seconds().loc[data_mean_v.idxmax()]
    p0 = [signal_max_initial, k_initial, t0_initial]

    if function in [_gompertz_drift, _sigmoid_drift]:
        p0 = p0 + [k_initial / 2, 2 * 60 * 360]

    params, _ = scipy.optimize.curve_fit(function, time, data[data_column], p0=p0)
    return params


if __name__ == "__main__":
    pr.plot_setup()
    data, ctrls = load_files()
    data = data[data.Read == 'GFP-G70-485,528']
    # cur_well_data = data[data.Well == "B14"]
    cur_well_data = data[data.Name == "ot_mix15_mm_rt"]
    # pr.kinetic_analysis_per_well()

    fit_function = {"gompertz": _gompertz, "gompertz_drift": _gompertz_drift,
                    "sigmoid": pr._sigmoid, "sigmoid_drift": _sigmoid_drift,
                    "weibull": _weibull, "weibull_drift": _weibull_drift}
    fit_function_name = "sigmoid_drift"
    params = fit_to_logistic(cur_well_data, function=fit_function[fit_function_name])
    # plot_kinetics(cur_well_data, params, fit_function=fit_function[fit_function_name], title=fit_function_name)

    # results = pr.kinetic_analysis(cur_well_data, group_by=["Name", "Read"],
    #                               fit_function_name=fit_function_name)
    # results["gain"] = results["Read"]
    # fit_function_name = "sigmoid"
    fit_function_name = "sigmoid_drift"
    results = pr.kinetic_analysis(cur_well_data, group_by=["Name", "Read"], fit_function_name=fit_function_name)
    pr.plot_kinetics(data, fit_function_name=fit_function_name)  #, kinetics=results)
    plt.show()

# Fit to kinetics data w/ pr

# Improve fits

# drift and background
# vectorize data fits
# add confidence intervals for fit curves
# Allow toggling b/w raw and smoothed averages
