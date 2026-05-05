import argparse
import ast
import json
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc

REQUIRED_BOUNDS_COLUMNS = {"Component", "Lower Bound", "Upper Bound"}


def _format_reagent_column_name(reagent, unit=None):
    reagent_name = str(reagent).strip()
    if not reagent_name:
        raise ValueError("Reagent name cannot be empty.")

    if unit is None or pd.isna(unit):
        return reagent_name

    unit_name = str(unit).strip()
    if not unit_name:
        return reagent_name
    if reagent_name.endswith(f" {unit_name}"):
        return reagent_name
    return f"{reagent_name} {unit_name}"


def load_bounds(df):
    """Load DOE bounds from a dataframe."""
    component_col = "Component"
    unit_col = _pick_column(df.columns, ("unit", "units"))
    names = df[component_col].astype(str).str.strip().tolist()
    if unit_col is not None:
        names = [
            _format_reagent_column_name(reagent=name, unit=unit)
            for name, unit in zip(names, df[unit_col])
        ]
    l_bounds = df["Lower Bound"].to_numpy(dtype=float)
    u_bounds = df["Upper Bound"].to_numpy(dtype=float)
    return names, l_bounds, u_bounds


def generate_doe(n_samples, l_bounds, u_bounds, names, log_space=False, seed=42):
    # LHS generator
    n_dims = len(l_bounds)
    if log_space:
        l_bounds = np.log(l_bounds)
        u_bounds = np.log(u_bounds)

    # LHS in [0, 1]^d
    sampler = qmc.LatinHypercube(d=n_dims, seed=seed)
    X_unit = sampler.random(n=n_samples)

    # Scale in bounds
    X = qmc.scale(X_unit, l_bounds=l_bounds, u_bounds=u_bounds)
    if log_space:  # Back-transform to original units
        X = np.exp(X)

    samples = pd.DataFrame(X, columns=names)
    return samples


def generate_ratio_sweep(
    ratios,
    n_total_concs,
    lower_bound,
    upper_bound,
    col_a="vioC uM",
    col_b="vioD uM",
    log_total=False,
):
    """Generate a 2D ratio × total-concentration sweep for two components.

    For each ratio r = col_a / col_b, the feasible total-concentration range
    [T_min, T_max] is computed from the per-enzyme bounds and sampled at
    n_total_concs evenly-spaced (or log-spaced) points.

    Args:
        ratios: list of col_a/col_b ratio values to sweep.
        n_total_concs: number of total-concentration levels per ratio.
        lower_bound: per-enzyme lower bound (same units as upper_bound).
        upper_bound: per-enzyme upper bound.
        col_a: output column name for the numerator component.
        col_b: output column name for the denominator component.
        log_total: if True, sample total concentrations in log space.

    Returns:
        DataFrame with columns [col_a, col_b], index named "conditions".
    """
    rows = []
    for r in ratios:
        T_min = max(lower_bound * (1 + r) / r, lower_bound * (1 + r))
        T_max = min(upper_bound * (1 + r) / r, upper_bound * (1 + r))
        if T_min >= T_max:
            continue
        if log_total:
            totals = np.exp(np.linspace(np.log(T_min), np.log(T_max), n_total_concs))
        else:
            totals = np.linspace(T_min, T_max, n_total_concs)
        for T in totals:
            rows.append({col_a: r * T / (1 + r), col_b: T / (1 + r)})
    df = pd.DataFrame(rows)
    df.index.name = "conditions"
    return df


def generate_standards_dilution_curve(standards_spec):
    """Generate a dilution curve for one or more pure compounds.

    For each compound, samples n_points concentrations between min_conc and
    max_conc (linear or log space). Each row has exactly one nonzero compound
    column; all other compound columns are zero.

    Args:
        standards_spec: list of dicts, each with keys:
            compound (str), unit (str), min_conc (float), max_conc (float),
            n_points (int), log_space (bool, optional, default False).

    Returns:
        DataFrame with columns [Type, <compound_a> <unit_a>, ...],
        Type="standard" for all rows, index named "conditions".
    """
    compound_cols = [f"{s['compound']} {s['unit']}" for s in standards_spec]
    rows = []
    for spec in standards_spec:
        col = f"{spec['compound']} {spec['unit']}"
        log_space = bool(spec.get("log_space", False))
        n_points = int(spec["n_points"])
        min_conc = float(spec["min_conc"])
        max_conc = float(spec["max_conc"])
        if log_space:
            concs = np.exp(np.linspace(np.log(min_conc), np.log(max_conc), n_points))
        else:
            concs = np.linspace(min_conc, max_conc, n_points)
        for c in concs:
            row = {"Type": "standard"}
            for compound_col in compound_cols:
                row[compound_col] = c if compound_col == col else 0.0
            rows.append(row)

    df = pd.DataFrame(rows)
    df.index.name = "conditions"
    return df


def concat_conditions_csvs(paths, default_type="doe"):
    """Concatenate multiple conditions CSVs into one combined CSV.

    Files without a ``Type`` column get
    ``Type = default_type``. The ``conditions`` index is
    reassigned 0..N-1 across all rows. Columns missing in some files are
    filled with 0.0.

    Args:
        paths: iterable of Path or str pointing to conditions CSVs.
        default_type: value written into ``Type`` for any file
            that does not already have that column.

    Returns:
        DataFrame with ``conditions`` as the index name, ``Type`` as
        the first column, followed by all other columns in order of first
        appearance.
    """
    dfs = []
    for path in paths:
        df = pd.read_csv(path)
        if "conditions" in df.columns:
            df = df.drop(columns=["conditions"])
        if "Type" not in df.columns:
            df.insert(0, "Type", default_type)
        dfs.append(df)

    combined = pd.concat(dfs, axis=0, ignore_index=True).fillna(0.0)

    cols = list(combined.columns)
    cols.insert(0, cols.pop(cols.index("Type")))
    combined = combined[cols]
    combined.index.name = "conditions"
    return combined


def visualize_doe(samples, names, log_scale=False):
    import matplotlib.pyplot as plt

    x = samples[names[0]]
    y = samples[names[1]]
    z = samples[names[2]]

    if log_scale:
        x = np.log10(x)
        y = np.log10(y)
        z = np.log10(z)

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(x, y, z, c=z, cmap="viridis", s=20, alpha=0.8)

    ax.set_xlabel(f'{"log (" + names[0] + ")" if log_scale else names[0]}')
    ax.set_ylabel(f'{"log (" + names[1] + ")" if log_scale else names[1]}')
    ax.set_zlabel(f'{"log (" + names[2] + ")" if log_scale else names[2]}')

    fig.colorbar(sc, label=f'{"log (" + names[2] + ")" if log_scale else names[2]}')
    plt.show()


def calc_reagent_vols(
    samples,
    final_volume_ul,
    reagents,
    used_volume,
    reagent_appendix="_conc_to_add",
):
    reagent_names = [c.replace("_conc_to_add", "") for c in samples.columns]
    reagent_names = list(dict.fromkeys(reagent_names))

    for reagent in reagent_names:
        stock = reagents.loc[reagent, "stock_conc"]  # scalar

        # Create stock_conc column for this reagent
        samples[f"{reagent}_stock_conc"] = stock

        # Compute volume to add (assuming same concentration units)
        # Handle stock = 0 or NaN safely
        assert pd.notna(stock) and stock > 0, "stock not correct"
        samples[f"{reagent}_vol_to_add"] = (
            samples[reagent + reagent_appendix] * final_volume_ul / stock
        )

    # check that all volumes are within the available dead volume
    all_added_reagents = [f"{reagent}_vol_to_add" for reagent in reagent_names]
    samples["water_vol_to_add"] = (
        final_volume_ul - used_volume - samples[all_added_reagents].sum(axis=1)
    )

    mask_bad = samples["water_vol_to_add"] < 0
    bad_idx = samples.index[mask_bad]
    if mask_bad.any():
        print(
            "WARNING: Lower sweep reagent concentrations -- not enough dead "
            "volume for some reaction conditions. "
            f"Offending sample indices: {list(bad_idx)}. Dropping these "
            "samples."
        )
        samples.drop(bad_idx, inplace=True)
    return samples


def calc_titration_concs(samples, final_rxn):
    for reagent in samples.columns:
        samples[f"{reagent}_conc_to_add"] = samples[reagent] - final_rxn.loc[reagent, "final_conc"]
    return samples


def _validate_bounds_df(bounds_df):
    missing = REQUIRED_BOUNDS_COLUMNS - set(bounds_df.columns)
    if missing:
        missing_display = ", ".join(sorted(missing))
        raise ValueError(
            "Bounds CSV is missing required column(s): "
            f"{missing_display}."
        )

    if bounds_df.empty:
        raise ValueError("Bounds CSV is empty.")

    if (bounds_df["Upper Bound"] < bounds_df["Lower Bound"]).any():
        raise ValueError(
            "Each 'Upper Bound' value must be greater than or equal to "
            "'Lower Bound'."
        )


def _normalize_column_name(name):
    return str(name).strip().lower().replace(" ", "_")


def _pick_column(columns, candidates):
    normalized = {_normalize_column_name(col): col for col in columns}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    return None


def _coerce_fixed_reagents_mapping(raw_mapping):
    if not isinstance(raw_mapping, dict):
        raise ValueError("Fixed reagents input must be a dictionary.")

    fixed_reagents = {}
    for key, value in raw_mapping.items():
        reagent = str(key).strip()
        if not reagent:
            raise ValueError("Fixed reagent names cannot be empty.")
        try:
            fixed_reagents[reagent] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Fixed reagent `{reagent}` has non-numeric value `{value}`."
            ) from exc

    if not fixed_reagents:
        raise ValueError("Fixed reagents dictionary is empty.")

    return fixed_reagents


def _parse_fixed_reagents_literal(raw_value):
    parsers = (json.loads, ast.literal_eval)
    parsed = None
    for parser in parsers:
        try:
            parsed = parser(raw_value)
            break
        except Exception:
            continue

    if parsed is None:
        raise ValueError(
            "Could not parse fixed reagents dictionary. Provide valid JSON/Python "
            "dict syntax or a CSV file path."
        )

    return _coerce_fixed_reagents_mapping(parsed)


def _load_fixed_reagents_csv(csv_path):
    df = pd.read_csv(csv_path)
    if df.empty:
        raise ValueError(f"Fixed reagents CSV is empty: {csv_path}")

    reagent_col = _pick_column(df.columns, ("reagent", "component", "name"))
    if reagent_col is not None:
        unit_col = _pick_column(df.columns, ("unit", "units"))
        value_col = _pick_column(
            df.columns,
            (
                "final_conc",
                "fixed_conc",
                "value",
                "concentration",
                "conc",
            ),
        )

        if value_col is None and len(df.columns) == 2:
            value_col = next(col for col in df.columns if col != reagent_col)

        if value_col is None:
            raise ValueError(
                "Could not identify value column in fixed reagents CSV. "
                "Include one of: final_conc, fixed_conc, value, concentration, "
                "or conc."
            )

        required_cols = [reagent_col, value_col]
        if unit_col is not None:
            required_cols.append(unit_col)
        fixed_df = df[required_cols].copy()
        fixed_df = fixed_df.dropna(subset=[reagent_col, value_col])
        fixed_df[reagent_col] = fixed_df[reagent_col].astype(str).str.strip()
        fixed_df["formatted_reagent"] = fixed_df.apply(
            lambda row: _format_reagent_column_name(
                reagent=row[reagent_col],
                unit=row[unit_col] if unit_col is not None else None,
            ),
            axis=1,
        )

        if fixed_df.empty:
            raise ValueError(
                f"Fixed reagents CSV has no valid reagent/value rows: {csv_path}"
            )

        duplicates = fixed_df["formatted_reagent"][
            fixed_df["formatted_reagent"].duplicated()
        ]
        if not duplicates.empty:
            duplicate_names = sorted(set(duplicates.tolist()))
            raise ValueError(
                "Duplicate reagent names found in fixed reagents CSV: "
                f"{duplicate_names}"
            )

        values = pd.to_numeric(fixed_df[value_col], errors="coerce")
        if values.isna().any():
            bad_rows = fixed_df.index[values.isna()].tolist()
            raise ValueError(
                "Non-numeric fixed reagent values found at rows "
                f"{bad_rows} in {csv_path}."
            )

        fixed_df[value_col] = values.astype(float)
        return dict(zip(fixed_df["formatted_reagent"], fixed_df[value_col]))

    if len(df) != 1:
        raise ValueError(
            "Fixed reagents CSV without a reagent column must have exactly one "
            "row (wide format with one column per reagent)."
        )

    row = df.iloc[0]
    fixed_reagents = {}
    for column, value in row.items():
        if pd.isna(value):
            continue
        reagent = _format_reagent_column_name(reagent=column, unit=None)
        if not reagent:
            continue
        try:
            fixed_reagents[reagent] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "Non-numeric value in fixed reagents CSV for column "
                f"`{column}`: `{value}`."
            ) from exc

    if not fixed_reagents:
        raise ValueError(
            f"Fixed reagents CSV has no usable reagent/value pairs: {csv_path}"
        )
    return fixed_reagents


def _resolve_fixed_reagents_path(raw_value, config_parent):
    raw_path = Path(raw_value).expanduser()
    if raw_path.is_absolute():
        return raw_path.resolve()

    cwd_candidate = (Path.cwd() / raw_path).resolve()
    if cwd_candidate.exists():
        return cwd_candidate
    return (config_parent / raw_path).resolve()


def _load_fixed_reagents(raw_value, config_parent):
    if raw_value is None:
        return {}

    candidate_path = _resolve_fixed_reagents_path(raw_value, config_parent)
    if candidate_path.exists():
        if not candidate_path.is_file():
            raise ValueError(
                f"Fixed reagents path must be a file: {candidate_path}"
            )
        return _load_fixed_reagents_csv(candidate_path)

    return _parse_fixed_reagents_literal(raw_value)


def _append_fixed_reagents(samples, fixed_reagents):
    if not fixed_reagents:
        return samples
    overlapping = sorted(set(samples.columns).intersection(fixed_reagents))
    if overlapping:
        raise ValueError(
            "Fixed reagent names overlap sweep columns: "
            f"{', '.join(overlapping)}"
        )
    for reagent, value in fixed_reagents.items():
        samples[reagent] = value
    return samples


def _apply_replicates(samples, replicates):
    if replicates < 1:
        raise ValueError("replicates must be a positive integer.")
    if replicates > 1:
        n_unique = len(samples)
        samples = samples.loc[samples.index.repeat(replicates)].reset_index(drop=True)
        samples.insert(0, "replicate", list(range(replicates)) * n_unique)
    samples.index = range(len(samples))
    samples.index.name = "conditions"
    return samples


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description="Generate DOE samples (Latin Hypercube, ratio sweep, or standards curve)."
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    # -- lhc subcommand (existing behaviour) --------------------------------
    lhc = subparsers.add_parser(
        "lhc",
        help=(
            "Generate Latin hypercube DOE samples from a bounds CSV and write "
            "the output CSV in the same directory as that bounds file."
        ),
    )
    lhc.add_argument(
        "bounds_csv",
        type=Path,
        help=(
            "Path to CSV with columns: Component, Lower Bound, Upper Bound. "
            "Optional unit/units column will be appended to output names as "
            "'<reagent> <unit>'."
        ),
    )
    lhc.add_argument(
        "n_samples",
        type=int,
        help="Number of DOE samples to generate.",
    )
    lhc.add_argument(
        "--fixed_reagents",
        nargs="?",
        default=None,
        help=(
            "Optional fixed reagents dictionary input. Provide either a dict "
            "string (for example '{\"hepes\": 40, \"atp\": 2}') or a CSV file "
            "path (for example 'fixed_final_rxn_concs.csv'). CSV can include "
            "optional unit/units column to produce '<reagent> <unit>' names."
        ),
    )
    lhc.add_argument(
        "--output-name",
        default=None,
        help=(
            "Optional output filename. If omitted, defaults to "
            "<bounds-stem>_doe_n<n_samples>.csv."
        ),
    )
    lhc.add_argument(
        "--log-space",
        action="store_true",
        help="Generate DOE in log space before back-transforming.",
    )
    lhc.add_argument(
        "--replicates",
        type=int,
        default=1,
        help="Number of replicates per DOE sample (default: 1, no replication).",
    )
    lhc.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for Latin hypercube sampling.",
    )

    # -- ratio-sweep subcommand ---------------------------------------------
    ratio = subparsers.add_parser(
        "ratio-sweep",
        help=(
            "Generate a 2D ratio × total-concentration sweep for two components. "
            "All parameters are read from a TOML config file."
        ),
    )
    ratio.add_argument(
        "config_toml",
        type=Path,
        help=(
            "Path to a TOML file defining the sweep. Required keys: "
            "col_a, col_b, lower_bound, upper_bound, ratios, n_total_concs. "
            "Optional keys: log_total (bool), fixed_reagents (path or dict), "
            "replicates (int), output_name (str)."
        ),
    )

    # -- standards-curve subcommand -----------------------------------------
    standards = subparsers.add_parser(
        "standards-curve",
        help=(
            "Generate dilution curves for one or more pure compounds. "
            "All parameters are read from a TOML config file."
        ),
    )
    standards.add_argument(
        "config_toml",
        type=Path,
        help=(
            "Path to a TOML file with [[standards]] entries. Required keys per "
            "entry: compound, unit, min_conc, max_conc, n_points. "
            "Optional keys: log_space (bool), output_name (str)."
        ),
    )

    # -- concat subcommand --------------------------------------------------
    concat = subparsers.add_parser(
        "concat",
        help=(
            "Concatenate multiple conditions CSVs into one. "
            "Files without a Type column are labelled with --default-type. "
            "Conditions are re-indexed 0..N-1 and missing columns filled with 0."
        ),
    )
    concat.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Two or more conditions CSV files to concatenate, in order.",
    )
    concat.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output CSV path. Defaults to combined_n<N>.csv written next to "
            "the first input file."
        ),
    )
    concat.add_argument(
        "--default-type",
        default="doe",
        help=(
            "Type value assigned to files that lack a Type column "
            "(default: doe)."
        ),
    )

    return parser.parse_args()


def _main_lhc(args):
    bounds_path = args.bounds_csv.expanduser().resolve()
    if not bounds_path.exists():
        raise FileNotFoundError(f"Bounds CSV not found: {bounds_path}")
    if not bounds_path.is_file():
        raise ValueError(f"Bounds path must be a file: {bounds_path}")

    if args.n_samples <= 0:
        raise ValueError("n_samples must be a positive integer.")

    bounds_df = pd.read_csv(bounds_path)
    _validate_bounds_df(bounds_df)

    names, l_bounds, u_bounds = load_bounds(bounds_df)
    samples = generate_doe(
        n_samples=args.n_samples,
        l_bounds=l_bounds,
        u_bounds=u_bounds,
        names=names,
        log_space=args.log_space,
        seed=args.seed,
    )

    fixed_reagents = _load_fixed_reagents(
        raw_value=args.fixed_reagents,
        config_parent=bounds_path.parent,
    )
    samples = _append_fixed_reagents(samples, fixed_reagents)
    samples = _apply_replicates(samples, args.replicates)

    output_name = args.output_name or (
        f"{bounds_path.stem}_doe_n{args.n_samples}.csv"
    )
    output_path = bounds_path.parent / output_name
    samples.to_csv(output_path, index=True)

    print(f"Wrote DOE samples to: {output_path}")
    if fixed_reagents:
        reagents_display = ", ".join(fixed_reagents.keys())
        print(f"Added fixed reagents: {reagents_display}")


def _main_ratio_sweep(args):
    config_path = args.config_toml.expanduser().resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Config TOML not found: {config_path}")

    with open(config_path, "rb") as f:
        cfg = tomllib.load(f)

    required_keys = {"col_a", "col_b", "lower_bound", "upper_bound", "ratios", "n_total_concs"}
    missing = required_keys - cfg.keys()
    if missing:
        raise ValueError(f"ratio_sweep_params.toml is missing required keys: {sorted(missing)}")

    col_a = cfg["col_a"]
    col_b = cfg["col_b"]
    lower_bound = float(cfg["lower_bound"])
    upper_bound = float(cfg["upper_bound"])
    ratios = [float(r) for r in cfg["ratios"]]
    n_total_concs = int(cfg["n_total_concs"])
    log_total = bool(cfg.get("log_total", False))
    replicates = int(cfg.get("replicates", 1))

    samples = generate_ratio_sweep(
        ratios=ratios,
        n_total_concs=n_total_concs,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        col_a=col_a,
        col_b=col_b,
        log_total=log_total,
    )

    fixed_reagents = _load_fixed_reagents(
        raw_value=cfg.get("fixed_reagents"),
        config_parent=config_path.parent,
    )
    samples = _append_fixed_reagents(samples, fixed_reagents)
    samples = _apply_replicates(samples, replicates)

    n_conditions = len(samples) // replicates if replicates > 1 else len(samples)
    default_output_name = f"{col_a.split()[0]}_{col_b.split()[0]}_ratio_sweep_n{n_conditions}.csv"
    output_name = cfg.get("output_name", default_output_name)
    output_path = config_path.parent / output_name
    samples.to_csv(output_path, index=True)

    print(f"Wrote ratio sweep ({len(ratios)} ratios × {n_total_concs} totals = {n_conditions} conditions) to: {output_path}")
    if fixed_reagents:
        reagents_display = ", ".join(fixed_reagents.keys())
        print(f"Added fixed reagents: {reagents_display}")


def _main_standards_curve(args):
    config_path = args.config_toml.expanduser().resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Config TOML not found: {config_path}")

    with open(config_path, "rb") as f:
        cfg = tomllib.load(f)

    if not cfg.get("standards"):
        raise ValueError(
            "standards_params.toml must contain at least one [[standards]] entry."
        )

    standards_spec = cfg["standards"]
    required_keys = {"compound", "unit", "min_conc", "max_conc", "n_points"}
    for i, spec in enumerate(standards_spec):
        missing = required_keys - spec.keys()
        if missing:
            raise ValueError(
                f"[[standards]] entry {i} is missing required keys: {sorted(missing)}"
            )

    # TOML quirk: keys placed after the last [[standards]] block land inside
    # that entry rather than at the top level. Check both places.
    _last = standards_spec[-1]
    replicates = int(cfg.get("replicates") or _last.get("replicates", 1))

    df = generate_standards_dilution_curve(standards_spec)
    df = _apply_replicates(df, replicates)

    n_unique = len(df) // replicates if replicates > 1 else len(df)
    output_name = cfg.get("output_name") or _last.get("output_name") or f"standards_n{n_unique}.csv"
    output_path = config_path.parent / output_name
    df.to_csv(output_path, index=True)

    replicate_note = f" × {replicates} replicates = {len(df)} rows" if replicates > 1 else ""
    print(f"Wrote {n_unique} standard conditions{replicate_note} to: {output_path}")
    for spec in standards_spec:
        scale = "log" if spec.get("log_space") else "linear"
        print(
            f"  {spec['compound']} {spec['unit']}: "
            f"{spec['n_points']} points [{spec['min_conc']}, {spec['max_conc']}] ({scale})"
        )


def _main_concat(args):
    paths = [p.expanduser().resolve() for p in args.inputs]
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(f"Input CSV not found: {p}")

    if len(paths) < 2:
        raise ValueError("concat requires at least two input files.")

    combined = concat_conditions_csvs(paths, default_type=args.default_type)

    if args.output is not None:
        output_path = args.output.expanduser().resolve()
    else:
        output_path = paths[0].parent / f"combined_n{len(combined)}.csv"

    combined.to_csv(output_path, index=True)

    counts = combined["Type"].value_counts().to_dict() if "Type" in combined.columns else {}
    counts_str = ", ".join(f"{v} {k}" for k, v in counts.items())
    print(f"Wrote {len(combined)} rows ({counts_str}) to: {output_path}")


def main():
    args = _parse_args()
    if args.mode == "lhc":
        _main_lhc(args)
    elif args.mode == "ratio-sweep":
        _main_ratio_sweep(args)
    elif args.mode == "standards-curve":
        _main_standards_curve(args)
    elif args.mode == "concat":
        _main_concat(args)


if __name__ == "__main__":
    main()
