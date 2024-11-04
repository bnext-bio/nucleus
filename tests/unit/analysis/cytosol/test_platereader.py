import pytest

PLATEMAP_CSV_DATA_PATH = "tests/test_data/platemap.csv"
PLATEMAP_XLSX_DATA_PATH = "tests/test_data/platemap.xlsx"

CYTATION_KINETIC_DATA_PATH = (
    "tests/test_data/cytation_fluorescence_kinetics.txt"
)
ENVISION_KINETIC_DATA_PATH = (
    "tests/test_data/envision_fluorescence_kinetics.csv"
)

DNA_SWEEP_DATA_PATH = "tests/test_data/cytation_dna_sweep.txt"
DNA_SWEEP_PLATEMAP_PATH = "tests/test_data/platemap.csv"


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


def test_read_cytation_temperature_column():
    from cdk.analysis.cytosol.platereader import read_cytation

    data = read_cytation(
        DNA_SWEEP_DATA_PATH
    )  # DNA sweep data has a weird temp degree character
    assert "Temperature (C)" in data.columns


@pytest.fixture(params=["envision", "cytation"])
def platereader_data(request):
    from cdk.analysis.cytosol.platereader import read_envision, read_cytation

    if request.param == "envision":
        return read_envision(ENVISION_KINETIC_DATA_PATH)
    elif request.param == "cytation":
        return read_cytation(CYTATION_KINETIC_DATA_PATH)


def test_platereader_read_has_mandatory_columns(platereader_data):
    data = platereader_data
    required_columns = {
        "Well",
        "Row",
        "Column",
        "Time",
        "Seconds",
        "Temperature (C)",
        "Read",
        "Data",
    }
    assert required_columns.issubset(
        set(data.columns)
    ), f"Missing required columns. Expected {required_columns} to be a subset of {set(data.columns)}"


def test_platereader_read_mandatory_columns_are_ordered(platereader_data):
    data = platereader_data
    required_columns = [
        "Well",
        "Row",
        "Column",
        "Time",
        "Seconds",
        "Temperature (C)",
        "Read",
        "Data",
    ]
    assert (
        required_columns == data.columns.values[: len(required_columns)]
    ).all(), f"Missing required columns. Expected {required_columns} to be a subset of {set(data.columns)}"


def test_platereader_time_is_timedelta(platereader_data):
    data = platereader_data
    print(data.dtypes)
    print(data["Time"].dtype)
    assert data["Time"].dtype == "timedelta64[ns]"


@pytest.fixture(
    params=[ENVISION_KINETIC_DATA_PATH, CYTATION_KINETIC_DATA_PATH]
)
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
    assert platemap.shape == (16, 5)


def test_read_platemap_normalizes_wells():
    import pandas as pd
    from cdk.analysis.cytosol.platereader import read_platemap

    platemap_check = pd.read_csv(PLATEMAP_CSV_DATA_PATH)
    assert (
        platemap_check["Well"].str.contains(":").any()
    ), "Test platemap does not contain wells with : (e.g., A:1); test is invalid."

    platemap = read_platemap(PLATEMAP_CSV_DATA_PATH)
    assert not platemap["Well"].str.contains(":").any()

@pytest.fixture()
def platemap_data(platemap_data_files):
    from cdk.analysis.cytosol.platereader import read_platemap

    return read_platemap(platemap_data_files)


def test_platemap_does_not_have_Float64_columns(platemap_data):
    import pandas as pd
    assert (pd.Float64Dtype() not in platemap_data.dtypes.values)
    # assert not any(col.dtype == pd.Float64Dtype() for col in data.columns)


def test_load_platereader_data_loads_platemap(
    platereader_data_file, platemap_data_files
):
    from cdk.analysis.cytosol.platereader import load_platereader_data

    data, platemap = load_platereader_data(
        platereader_data_file, platemap_data_files
    )

    # We will only have all of these columns if the platemap was successfully loaded.
    required_columns = {"Well", "Data", "DNA Template"}
    assert required_columns.issubset(set(data.columns))


def test_load_platereader_with_none_platemap(platereader_data_file):
    from cdk.analysis.cytosol.platereader import load_platereader_data

    data = load_platereader_data(platereader_data_file, None)
    assert data.shape[0] > 0


@pytest.fixture
def dna_sweep_data():
    from cdk.analysis.cytosol.platereader import load_platereader_data

    return load_platereader_data(
        DNA_SWEEP_DATA_PATH, DNA_SWEEP_PLATEMAP_PATH
    ).data


def test_plot_plate(dna_sweep_data):
    from cdk.analysis.cytosol.platereader import plot_plate

    g = plot_plate(dna_sweep_data)
    assert g.axes.shape[0] > 0 and g.axes.shape[1] > 0
