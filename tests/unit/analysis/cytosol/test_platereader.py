import pytest

PLATEMAP_CSV_DATA_PATH = "tests/test_data/platemap.csv"
PLATEMAP_XLSX_DATA_PATH = "tests/test_data/platemap.xlsx"

CYTATION_KINETIC_DATA_PATH = "tests/test_data/cytation_fluorescence_kinetics.txt"
ENVISION_KINETIC_DATA_PATH = "tests/test_data/envision_fluorescence_kinetics.csv"


def test_read_envision():
    import pandas as pd
    from cdk.analysis.cytosol.platereader import read_envision

    data = read_envision(ENVISION_KINETIC_DATA_PATH)
    print(f"Columns: {data.columns}")

    assert isinstance(data, pd.DataFrame)


def test_read_cytation():
    import pandas as pd
    from cdk.analysis.cytosol.platereader import read_cytation

    data = read_cytation(CYTATION_KINETIC_DATA_PATH)
    print(f"Columns: {data.columns}")

    assert isinstance(data, pd.DataFrame)


@pytest.fixture(params=["envision", "cytation"])
def platereader_data(request):
    from cdk.analysis.cytosol.platereader import read_envision, read_cytation

    if request.param == "envision":
        return read_envision(ENVISION_KINETIC_DATA_PATH)
    elif request.param == "cytation":
        return read_cytation(CYTATION_KINETIC_DATA_PATH)


def test_platereader_read_has_mandatory_columns(platereader_data):
    data = platereader_data
    required_columns = {"Well", "Row", "Column", "Time", "Seconds", "Temperature (C)", "Read", "Data"}
    assert required_columns.issubset(
        set(data.columns)
    ), f"Missing required columns. Expected {required_columns} to be a subset of {set(data.columns)}"


def test_platereader_read_mandatory_columns_are_ordered(platereader_data):
    data = platereader_data
    required_columns = ["Well", "Row", "Column", "Time", "Seconds", "Temperature (C)", "Read", "Data"]
    assert (
        required_columns == data.columns.values[: len(required_columns)]
    ).all(), f"Missing required columns. Expected {required_columns} to be a subset of {set(data.columns)}"


@pytest.fixture(params=[ENVISION_KINETIC_DATA_PATH, CYTATION_KINETIC_DATA_PATH])
def platereader_data_file(request):
    return request.param


def test_load_platereader_data(platereader_data_file):
    import pandas as pd
    from cdk.analysis.cytosol.platereader import load_platereader_data

    data = load_platereader_data(platereader_data_file)

    assert isinstance(data, pd.DataFrame)


@pytest.fixture(params=[PLATEMAP_CSV_DATA_PATH, PLATEMAP_XLSX_DATA_PATH])
def platemap_data_files(request):
    return request.param


def test_read_platemap(platemap_data_files):
    from cdk.analysis.cytosol.platereader import read_platemap

    platemap = read_platemap(platemap_data_files)
    print(platemap)
    assert platemap.shape == (16, 3)
