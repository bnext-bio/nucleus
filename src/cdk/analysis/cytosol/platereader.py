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
import glob
import datetime

import numpy as np
import pandas as pd

def read_glomax(data_dir: str) -> pd.DataFrame:
    all_df = list()

    # glob over .csv files in dfpath; append to data; concatenate into one DataFrame
    for csv in glob.glob(f"{data_dir}/*.csv"):
        df = pd.read_csv(csv)
        df["File"] = os.path.basename(csv)
        df["Row"] = df["WellPosition"].str.split(":").str[0]
        df["Column"] = df["WellPosition"].str.split(":").str[1].astype(int) 
        df["Time"] = pd.to_datetime(data["File"].str.replace(r'.* ([0-9.]+ [0-9_]+).*', r'\1', regex=True), format="%Y.%m.%d %H_%M_%S")
        df["WellTime"] = pd.to_timedelta(data["Timestamp(ms)"], "us")
    
        data.append(df) 
    
    data = pd.concat(data, ignore_index=True)

    # label different wavelengths
    channel_map = dict(zip(data["ID"].unique(), ["A600", "Blue", "Green", "Red"]))
    data["Channel"] = data["ID"].map(channel_map)

    palette = dict(zip(dict.fromkeys(channel_map.values()), ["brown", "limegreen", "red", "firebrick"]))

    # massage Time
    data["TimeDelta"] = data["Time"] - data["Time"].min()
    data["TimeDeltaPretty"] = data["TimeDelta"].map(lambda x: "{:02d}:00".format(x.components.hours)) #{:02d}".format(x.components.hours, x.components.minutes))

    # Get a generic data column
    data["Data"] = data["CalculatedFlux"]
    data["Data"].fillna(data[data["CalculatedFlux"].isna()]["OpticalDensity"], inplace=True)
    
    # Label replicates
    data["Replicate"] = data["File"].map(lambda x: re.sub(r'.* OUT ([0-9]+).csv', r'\1', x))
    data.sort_values(by=["TimeDelta", "Row", "Column"], inplace=True)

    return data

def read_cytation(data_file: str, sep="\t") -> pd.DataFrame:
    # read data file as long string
    data = ""
    with open(data_file, "r", encoding="latin1") as file:
        data = file.read()

    # extract indices for Proc Details, Layout
    procidx = re.search(r"Procedure Details", data)
    layoutidx = re.search(r"Layout", data)
    readidx = re.search(r"^(Read\s)?\d+,\d+", data, re.MULTILINE)

    # get header DataFrame
    header = data[:procidx.start()]
    header = pd.read_csv(
        io.StringIO(header), 
        delimiter=sep, 
        header=0, names=["key", "value"]
    )

    # get procedure DataFrame
    procedure = data[procidx.end():layoutidx.start()]  
    procedure = pd.read_csv(io.StringIO(procedure), skipinitialspace=True, names=range(4))
    procedure = procedure.replace(np.nan, "")

    # get Cytation plate map from data_file as DataFrame
    layout = data[layoutidx.end():readidx.start()]
    layout = pd.read_csv(io.StringIO(layout), index_col=False)
    layout = layout.set_index(layout.columns[0])
    layout.index.name = "Row"

    # iterate over data string to find individual reads
    reads = dict()
    
    sep = r"(?:Read\s\d+:)?(?:\s\d{3},\d{3}(?:\[\d\])?)?" + sep
        
    for readidx in re.finditer(r"^(Read\s)?\d+,\d+.*\n", data, re.MULTILINE):
        # for each iteration, extract string from start idx to end icx
        read = data[readidx.end():]
        read = read[:re.search(r"(^(Read\s)?\d+,\d+|^Blank Read\s\d|Results|\Z)", read[1:], re.MULTILINE).start()]
        read = pd.read_csv(io.StringIO(read), sep=sep, engine="python").convert_dtypes()
        reads[data[readidx.start():readidx.end()].strip()] = read

    # create a DataFrame for each read and process, then concatenate into a large DataFrame
    # NOTE: JC 2024-05-21 - turns out, len(list(reads.items())) = 1 (one big mono table)
    read_dataframes = list()
    for name, r in reads.items():
        # filter out Cytation calculated kinetic parameters, which are cool, but don't want rn
        r = r[r.Time.str.contains(r"\d:\d{2}:\d{2}", regex=True)]
        
        # extract meaningful parameters from really big string
        r = r.melt(id_vars=["Time", "T°"], var_name="Well", value_name="Data")
        r['Row'] = r['Well'].str.extract(r"([A-Z]+)")
        r['Column'] = r['Well'].str.extract(r"(\d+)").astype(int)
        r['Temperature (C)'] = r['T°'].str.extract(r"(\d+)").astype(float)
        r['Data'] = r['Data'].replace("OVRFLW", np.inf)
        r['Data'] = r['Data'].astype(float)
        r['Read'] = name
        r['Ex'] = r['Read'].str.extract(r"(\d+),\d+").astype(int)
        r['Em'] = r['Read'].str.extract(r"\d+,(\d+)").astype(int)
        read_dataframes.append(r)

    data = pd.concat(read_dataframes)

    # add time column to data DataFrame
    data['Time'] = pd.to_timedelta(data['Time']).astype('timedelta64[s]')
    data['Seconds'] = data['Time'].map(lambda x: x.total_seconds())

    return data[['Well', 'Row', 'Column', 'Time', 'Seconds', 'Temperature (C)', 'Read', 'Data']]

def read_envision(data_file: str) -> pd.DataFrame:
    # load data
    data = pd.read_csv(data_file).convert_dtypes()

    # massage Row, Column, and Well information
    data["Row"] = data["Well ID"].apply(lambda s: s[0]).astype(pd.StringDtype())
    data["Column"] = data["Well ID"].apply(lambda s: str(int(s[1:])))
    data["Well"] = data.apply(lambda well: f"{well['Row']}:{well['Column']}", axis=1)

    data['Time'] = pd.to_timedelta(data['Time [hhh:mm:ss.sss]']).astype('timedelta64[s]')
    data['Seconds'] = data['Time'].map(lambda x: x.total_seconds())

    data['Temperature (C)'] = data['Temperature current[°C]']

    data['Read'] = data['Operation']

    data["Data"] = data["Result Channel 1"]

    data["Excitation (nm)"] = data["Exc WL[nm]"]
    data["Emission (nm)"] = data["Ems WL Channel 1[nm]"]
    data["Wavelength (nm)"] = data["Excitation (nm)"] + "," + data["Emission (nm)"]
        
    return data[['Well', 'Row', 'Column', 'Time', 'Seconds', 'Temperature (C)', 'Read', 'Data']]
