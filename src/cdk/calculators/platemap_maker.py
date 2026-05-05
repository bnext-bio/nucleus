import string
import pandas as pd
import random
import math

def generate_384_well_ids(order: str = "row", skip: bool = False,
                          randomize: bool = False,
                          rng: random.Random | None = None) -> list[str]:
    """
    Generate well IDs for a 384-well plate (rows A–P, cols 1–24).

    Parameters
    ----------
    order : {'row', 'column'}, default 'row'
        'row' -> A1, A2, ..., A24, B1, ...
        'column' -> A1, B1, ..., P1, A2, ...
    skip : bool, default False
        If True, keep every other well in the sequence.
    randomize : bool, default False
        If True, randomly shuffle the resulting list of wells.
    rng : random.Random or None, default None
        Optional RNG instance for reproducible shuffling. If None and
        shuffle is True, uses the global random module.
    """
    rows = list(string.ascii_uppercase[:16])   # A-P
    cols = list(range(1, 25))                  # 1-24

    wells: list[str] = []
    if order == "row":
        wells = [f"{r}{c}" for r in rows for c in cols]
    elif order == "column":
        wells = [f"{r}{c}" for c in cols for r in rows]
    else:
        raise ValueError("order must be 'row' or 'column'")

    if skip:
        wells = wells[::2]

    if randomize:
        if rng is None:
            random.shuffle(wells)
        else:
            rng.shuffle(wells)

    return wells

def add_well_ids_column(df: pd.DataFrame,
                        order: str = "row",
                        skip: bool = False,
                        randomize: bool = False,
                        rng: random.Random | None = None,
                        column_name: str = "Well") -> pd.DataFrame:
    """
    Add 384-well plate IDs as a column to a DataFrame.

    This function generates well identifiers for a 384-well plate
    (rows A–P, columns 1–24) in a specified traversal order and
    assigns them to the rows of the input DataFrame as a new column.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame to which well IDs will be assigned. The number
        of rows in the DataFrame must not exceed the number of available
        wells generated under the chosen options.
    order : {'row', 'column'}, default 'row'
        Ordering in which wells are traversed:
        - 'row': A1, A2, ..., A24, B1, B2, ..., P24
        - 'column': A1, B1, ..., P1, A2, B2, ..., P24
    skip : bool, default False
        If True, keep every other well in the generated sequence
        (i.e. use wells[::2]). This is useful when spacing samples
        to reduce edge or interaction effects.
    randomize : bool, default False
        If True, randomly shuffle the well IDs before assignment.
    rng : random.Random or None, default None
        Optional RNG for reproducible shuffling.
    column_name : str, default 'well_id'
        Name of the column to create in the returned DataFrame that
        will store the well identifiers.

    Returns
    -------
    pandas.DataFrame
        A copy of the input DataFrame with an additional column
        containing the assigned well IDs.

    Raises
    ------
    ValueError
        If the DataFrame has more rows than the number of wells
        available under the selected ordering and skip settings.

    Notes
    -----
    The generated well IDs follow the conventional 384-well plate
    layout with rows labeled A–P and columns labeled 1–24.
    """
    wells = generate_384_well_ids(order=order, skip=skip, randomize=randomize, rng=rng)
    if len(df) > len(wells):
        raise ValueError("DataFrame has more rows than available wells.")
    df = df.copy()
    df[column_name] = wells[:len(df)]
    return df


def _choose_centered_rectangle(
        n_samples: int,
        plate_rows: int = 16,
        plate_cols: int = 24,
) -> tuple[int, int]:
    """Pick a compact centered rectangle that can hold n_samples wells."""
    if n_samples <= 0:
        raise ValueError("n_samples must be > 0.")
    if n_samples > plate_rows * plate_cols:
        raise ValueError(
            f"n_samples ({n_samples}) exceeds plate capacity ({plate_rows * plate_cols})."
        )

    target_ratio = plate_rows / plate_cols
    best_dims: tuple[int, int] | None = None
    best_score: tuple[float, float, float] | None = None

    for rows_used in range(1, plate_rows + 1):
        cols_used = math.ceil(n_samples / rows_used)
        if cols_used > plate_cols:
            continue

        area = rows_used * cols_used
        ratio = rows_used / cols_used
        # Prefer smallest area first, then plate-like shape, then near-square.
        score = (area, abs(ratio - target_ratio), abs(rows_used - cols_used))
        if best_score is None or score < best_score:
            best_score = score
            best_dims = (rows_used, cols_used)

    if best_dims is None:
        raise ValueError(
            f"Could not place {n_samples} samples on a {plate_rows}x{plate_cols} plate."
        )
    return best_dims


def generate_centered_384_well_ids(
        n_samples: int,
        randomize: bool = True,
        rng: random.Random | None = None,
) -> list[str]:
    """Generate centered 384-well IDs from a compact middle rectangle.

    Wells are selected from a centered rectangle (rows A-P, columns 1-24),
    then optionally randomized. This keeps samples away from plate edges while
    preserving randomization for spatial bias reduction.
    """
    rows = list(string.ascii_uppercase[:16])  # A-P
    cols = list(range(1, 25))                 # 1-24

    rows_used, cols_used = _choose_centered_rectangle(
        n_samples=n_samples,
        plate_rows=len(rows),
        plate_cols=len(cols),
    )
    row_start = (len(rows) - rows_used) // 2
    col_start = (len(cols) - cols_used) // 2

    row_subset = rows[row_start: row_start + rows_used]
    col_subset = cols[col_start: col_start + cols_used]
    wells = [f"{r}{c}" for r in row_subset for c in col_subset]

    if randomize:
        if rng is None:
            random.shuffle(wells)
        else:
            rng.shuffle(wells)

    return wells[:n_samples]


def add_centered_well_ids_column(
        df: pd.DataFrame,
        randomize: bool = True,
        rng: random.Random | None = None,
        column_name: str = "Well",
) -> pd.DataFrame:
    """Add centered, randomized 384-well IDs to a DataFrame."""
    wells = generate_centered_384_well_ids(
        n_samples=len(df),
        randomize=randomize,
        rng=rng,
    )
    out = df.copy()
    out[column_name] = wells
    return out
