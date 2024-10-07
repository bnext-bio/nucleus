import pytest
import os.path

CYTATION_KINETIC_DATA_PATH = "tests/test_data/cytation_fluorescence_kinetics.txt"
ENVISION_KINETIC_DATA_PATH = "tests/test_data/envision_fluorescence_kinetics.csv"

@pytest.fixture(params=['envision', 'cytation'])
def platereader_data(request):
    from cdk.analysis.cytosol.platereader import read_envision, read_cytation
    if request.param == 'envision':
        return read_envision(ENVISION_KINETIC_DATA_PATH)
    elif request.param == 'cytation':
        return read_cytation(CYTATION_KINETIC_DATA_PATH)


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

def test_platereader_read_has_mandatory_columns(platereader_data):
    data = platereader_data
    required_columns = {"Well", "Row", "Column", "Time", "Seconds", "Temperature (C)", "Read", "Data"}
    assert required_columns.issubset(set(data.columns)), f"Missing required columns. Expected {required_columns} to be a subset of {set(data.columns)}"

def test_platereader_read_mandatory_columns_are_ordered(platereader_data):
    data = platereader_data
    required_columns = ["Well", "Row", "Column", "Time", "Seconds", "Temperature (C)", "Read", "Data"]
    assert (required_columns == data.columns.values[:len(required_columns)]).all(), f"Missing required columns. Expected {required_columns} to be a subset of {set(data.columns)}"
    # assert data.shape == (96, 18)
    # assert set(data.columns) == set([
    #     'Well ID', 'Result Channel 1', 'Exc WL[nm]', 'Ems WL Channel 1[nm]', 
    #     'Time [hhh:mm:ss.sss]', 'Repeat', 'Row', 'Column', 'Well', 'Data', 
    #     'Ex', 'Em', 'Wavelength', 'Time', 'Seconds', 'Background', 
    #     'BackgroundSubtracted'
    # ])

    # assert data['Data'].dtype == float
    # assert data['Time'].dtype == 'timedelta64[s]'
    # assert data['Seconds'].dtype == float