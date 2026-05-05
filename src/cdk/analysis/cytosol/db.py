import datetime
import logging
import pandas as pd
from pathlib import Path
import os
import re
from tempfile import TemporaryDirectory

import cdk
from cdk.analysis.cytosol import platereader as pr

# TODO:
#   - Properly track / handle the normalized data column so we can check and name it appropriately.

DB_DATA_DIR = "data/"
DB_SUMMARY_DIR = "summary/"
DB_COMPOSITE_DIR = "composite/"

DB_DATA_FILE_PREFIX = "db-cytosol-data"
DB_SUMMARY_FILE_PREFIX = "db-cytosol-summary"
DB_COMPOSITE_FILE_PREFIX = "db-cytosol"

CDK_VERSION_COL = "CDK Version"

REQUIRED_DATA_COLUMNS = [
    "Date",
    "Experiment",
    "Name",
    "Well",
    "Name",
    "Type",
    "data_normalized",
    CDK_VERSION_COL
]

REQUIRED_KINETICS_COLUMNS = [
    "Date",
    "Experiment",
    "Name",
    "Read",
    "Well",
    "Type"
]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class CytosolDB():
    def __init__(self, db_path: str = None):
        if db_path is not None:
            self._db_dir = Path(db_path)
        else:
            self._db_dir = Path(os.getenv("CDK_DB_CYTOSOL", "./db/"))

        log.info(f"Cytosol DB Initialized: database at {self._db_dir}")

    def load(self):
        log.info(f"Loading cytosol DB at {self._db_dir}")
        
        data_dir_path = self._db_dir / DB_DATA_DIR
        all_data = dict()
        
        for data_file in data_dir_path.glob(f"{DB_DATA_FILE_PREFIX}*.csv"):
            log.debug(f"Reading data file {data_file}")
            data = pd.read_csv(data_file)
            all_data[data_file] = data
            
        summary_dir_path = self._db_dir / DB_SUMMARY_DIR
        all_summary = dict()
        for summary_file in summary_dir_path.glob(f"{DB_SUMMARY_FILE_PREFIX}*.csv"):
            log.debug(f"Reading summary file {summary_file}")
            summary = pd.read_csv(summary_file)
            all_summary[summary_file] = summary
        
        log.info(f"Loaded {len(all_data)} data files and {len(all_summary)} summaries.")
        
        data = pd.concat(all_data.values()).reset_index(drop=True)
        data["Date"] = pd.to_datetime(data["Date"], format="mixed")
        data = data.sort_values(by=["Date", "Experiment"])
        
        summary = pd.concat(all_summary.values()).reset_index(drop=True)
        summary["Date"] = pd.to_datetime(summary["Date"], format="mixed")
        summary = summary.sort_values(by=["Date", "Experiment"])
        
        self._data = data
        self._summary = summary
        
        return data, summary

    def to_composite(self):
        composite_path = self._db_dir / DB_COMPOSITE_DIR
        composite_path.mkdir(parents=True, exist_ok=True)
        
        with TemporaryDirectory(dir=composite_path) as tmp:
            data_path, summary_path = self.to_csv(Path(tmp) / DB_COMPOSITE_FILE_PREFIX)
            data_path.replace(composite_path / f"{DB_COMPOSITE_FILE_PREFIX}-data.csv")
            summary_path.replace(composite_path / f"{DB_COMPOSITE_FILE_PREFIX}-summary.csv")
        
        log.info(f"Updated composite database at {composite_path}.")

    def to_csv(self, filename: str | Path, overwrite: bool = False):
        file_path = Path(filename)
        data_path = file_path.parent / f"{file_path.stem}-data.csv"
        summary_path = file_path.parent / f"{file_path.stem}-summary.csv"
        
        if not overwrite and (data_path.exists() or summary_path.exists()):
            log.warning(f"Output file exists, use overwrite=True to overwrite: {data_path}, {summary_path}.")
            return
        
        self._data.to_csv(data_path, index=False)
        self._summary.to_csv(summary_path, index=False)
        
        return data_path, summary_path
 
    def to_excel(self, filename: str | Path, overwrite: bool = False):
        file_path = Path(filename)
        data_path = file_path.parent / f"{file_path.stem}-data.xlsx"
        summary_path = file_path.parent / f"{file_path.stem}-summary.xlsx"
        
        if not overwrite and (data_path.exists() or summary_path.exists()):
            log.warning(f"Output file exists, use overwrite=True to overwrite: {data_path}, {summary_path}.")
            return
        
        self._data.to_excel(data_path)
        self._summary.to_excel(summary_path)
        
    def add_version(self, data: pd.DataFrame) -> pd.DataFrame:
        data[CDK_VERSION_COL] = cdk.__version__
        return data
    
    def validate_columns(self, df: pd.DataFrame, required_columns: list = REQUIRED_DATA_COLUMNS) -> bool:        
        for col in required_columns:
            if col not in df.columns:
                log.warning(f"Required column {col} not present in data (columns={", ".join(df.columns)}))")
                return False
            
        return True
        
    def export(self, data: pd.DataFrame, platemap: pd.DataFrame, kinetics: pd.DataFrame, force=False):
        data = data.rename(columns={"Experiment Name": "Experiment"}) 
        experiment = data["Experiment"].iloc[0]
        
        log.debug(f"Exporting experiment: {experiment}")
        exp_date = data["Clock Time"].min().strftime("%Y%m%d")
        exp_name = re.sub(r"\W+", "-", experiment.lower())
        exp_filename = f"{DB_DATA_FILE_PREFIX}-{exp_date}-{exp_name}.csv"
        log.debug(f"Experiment data file: {exp_filename}")

        data_path = self._db_dir / DB_DATA_DIR / exp_filename
        data_path.parent.mkdir(parents=True, exist_ok=True)
        log.debug(f"Output data path: {data_path}")

        summary_filename = f"{DB_SUMMARY_FILE_PREFIX}-{exp_date}-{exp_name}.csv"
        summary_path = self._db_dir / DB_SUMMARY_DIR / summary_filename
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        log.debug(f"Kinetics summary path: {summary_path}")
        
        if data_path.exists() or summary_path.exists():
            if force:
                log.info(f"Overwriting existing data file(s) {data_path} {summary_path}.")
            else:
                log.warning(f"Output file(s) already exist. Use `force=True` to overwrite ({data_path, summary_path}).")
                return
        
        data = self.add_version(data)
            
        if not self.validate_columns(data):
            log.warning("Data validation failed: will not export.")
            return
            
        data.to_csv(data_path, index=False)

        kinetics = kinetics.copy()
        kinetics.columns = [re.sub(r"_$", "", '_'.join(col)) for col in kinetics.columns]
        kinetics = kinetics.reset_index()
        
        # Remove common columns between kinetics and platemap,
        # except for Well which we'll merge on. This stops duplicates appearing
        # when those parts of the platemap are in the kinetics index because we
        # grouped on them.
        for col in platemap.columns:
            if col != "Well" and col in kinetics.columns:
                kinetics = kinetics.drop(columns=col)
        
        kinetics = platemap.merge(kinetics, how="right", on="Well")
        kinetics = kinetics.rename(columns={"Experiment Name": "Experiment"})
        kinetics = self.add_version(kinetics)
        
        if not self.validate_columns(kinetics, REQUIRED_KINETICS_COLUMNS):
            log.warning("Kinetics validation failed: will not export.")
            return

        kinetics.to_csv(summary_path, index=False)