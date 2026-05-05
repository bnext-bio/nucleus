# Discovery Plate Pipeline

End-to-end guide for designing experiment conditions with `doe.py` and computing
liquid-handling volumes with `discovery_plate_concentration_pipeline.py`.

All commands below assume you are running from the `cdk/` directory.

---

## Overview

```
bounds.csv / ratio_sweep_params.toml / standards_params.toml
          │
          ▼
       doe.py              ← generate conditions (lhc, ratio-sweep, standards-curve, concat)
          │
          ▼
 <design>_n<N>.csv         ← one row per condition, columns = "<reagent> <unit>"
          │
          ▼
discovery_plate_concentration_pipeline.py   ← compute master mix + titration volumes
          │
          ▼
 samples_titration.csv         ← per-well pipetting instructions
 samples_titration_labcraft.csv ← well-indexed for liquid handlers
 base_master_mix.csv           ← pooled base recipe
```

---

## Step 1: Generate Conditions with `doe.py`

`doe.py` has four subcommands: `lhc`, `ratio-sweep`, `standards-curve`, and `concat`.

### Latin Hypercube (`lhc`)

Randomly samples the N-dimensional parameter space defined by a bounds CSV.

**Bounds CSV format** (`bounds.csv`):

| Component | Lower Bound | Upper Bound | Unit |
|-----------|-------------|-------------|------|
| pmix | 1.2 | 2.7 | mg/ml |
| ribosome | 1.2 | 2.7 | uM |
| magnesium_acetate | 6 | 18 | mM |

`Unit` is optional; when present it is appended to the column name as `<reagent> <unit>` in the output.

**Command:**

```bash
PYTHONPATH=src python src/cdk/calculators/doe.py lhc \
  examples/discovery_plate_pipeline/bounds.csv \
  20 \
  --fixed_reagents examples/discovery_plate_pipeline/fixed_final_rxn_concs.csv
```

This writes `bounds_doe_n20.csv` next to `bounds.csv`.

**All `lhc` flags:**

| Flag | Default | Description |
|------|---------|-------------|
| `bounds_csv` | *(required)* | Path to bounds CSV |
| `n_samples` | *(required)* | Number of conditions to generate |
| `--fixed_reagents` | — | CSV file path or inline dict of reagents held constant across all conditions (e.g. `'{"hepes": 50}'`) |
| `--output-name` | `<bounds-stem>_doe_n<N>.csv` | Override output filename |
| `--log-space` | false | Sample in log space then back-transform (useful for orders-of-magnitude ranges) |
| `--replicates` | 1 | Repeat each condition N times (adds a `replicate` column) |
| `--seed` | 42 | Random seed |

---

### Ratio Sweep (`ratio-sweep`)

Generates a 2D grid of **ratio × total concentration** for two components.
For each ratio `r = col_a / col_b`, the feasible total-concentration range
within the per-enzyme bounds is computed and sampled at `n_total_concs` points.
This is more interpretable than LHC when the balance between two competing
components is the primary variable of interest.

All parameters are read from a TOML config file.

**`ratio_sweep_params.toml`:**

```toml
col_a = "vioC uM"           # numerator component (must match reagents.csv name)
col_b = "vioD uM"           # denominator component
lower_bound = 0.1           # per-enzyme lower bound (µM)
upper_bound = 3.0           # per-enzyme upper bound (µM)
ratios = [0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]   # col_a / col_b
n_total_concs = 5           # total-concentration levels per ratio
log_total = false           # if true, sample totals in log space
fixed_reagents = "fixed_final_rxn_concs.csv"        # optional, relative path
output_name = "bounds_vioC_vioD_ratio_sweep_n35.csv"
# replicates = 1            # optional
```

**Command:**

```bash
PYTHONPATH=src python src/cdk/calculators/doe.py ratio-sweep \
  path/to/experiment_setup/ratio_sweep_params.toml
```

Paths in `ratio_sweep_params.toml` (e.g. `fixed_reagents`) are resolved
relative to the TOML file itself, so the config is fully self-contained.

**Feasible range per ratio:**

For ratio `r` and per-enzyme bounds `[lo, hi]`:

```
T_min = max(lo*(1+r)/r,  lo*(1+r))
T_max = min(hi*(1+r)/r,  hi*(1+r))
```

---

### Standards Curve (`standards-curve`)

Generates a dilution series for one or more pure compounds. Each row has one
nonzero compound column; all others are zero. Used to co-plate calibration
standards alongside DOE conditions.

**`standards_params.toml`:**

```toml
[[standards]]
compound = "violacein"
unit = "uM"
min_conc = 6.5
max_conc = 100.0
n_points = 6
log_space = true

[[standards]]
compound = "deoxyviolacein"
unit = "uM"
min_conc = 6.5
max_conc = 100.0
n_points = 6
log_space = true

output_name = "standards_n12.csv"
# replicates = 1   # optional — place before or after [[standards]] blocks
```

> **Pipetting constraint:** `min_conc` must satisfy
> `min_conc ≥ min_pipetting_vol_ul × stock_conc / final_rxn_vol_ul`.
> With 125 µM stock, 20 µL final volume, and 1.0 µL minimum, the floor is 6.25 µM.

**Command:**

```bash
PYTHONPATH=src python src/cdk/calculators/doe.py standards-curve \
  path/to/experiment_setup/standards_params.toml
```

Output has `Type = "standard"` and one column per compound.
Add each compound to `reagents.csv` with its stock concentration before running the pipeline.

---

### Concat (`concat`)

Concatenates multiple conditions CSVs into one. Files without a `Type` column
get `Type = "doe"` (overrideable with `--default-type`). Conditions are
re-indexed 0..N-1 and missing columns filled with `0.0`.

```bash
PYTHONPATH=src python src/cdk/calculators/doe.py concat \
  bounds_vioC_vioD_ratio_sweep_n35.csv \
  standards_n12.csv \
  --output bounds_with_standards.csv
```

| Flag | Default | Description |
|------|---------|-------------|
| `inputs` | *(required)* | Two or more conditions CSV files, in order |
| `--output` | `combined_n<N>.csv` next to first input | Output CSV path |
| `--default-type` | `doe` | `Type` value for files that lack that column |

---

### `fixed_reagents` CSV format

Used by `lhc` and `ratio-sweep` to append constant columns to every condition:

```csv
reagent, final_conc, units
vioA, 0.1, uM
vioB, 1.0, uM
vioE, 1.0, uM
```

The `units` column is optional. When present, column names in the output are
formatted as `<reagent> <unit>` (e.g. `vioA uM`).

---

## Step 2: Run the Concentration Pipeline

Once you have a conditions CSV, the pipeline computes:
- A **base master mix** of reagents that are constant across all conditions
- **Per-well titration volumes** for reagents that vary
- For standard rows (`Type = "standard"`): simple `compound_vol + buffer_vol = final_rxn_vol_ul` — no master mix

**Command:**

```bash
PYTHONPATH=src python src/cdk/calculators/discovery_plate_concentration_pipeline.py \
  --experiment-dir path/to/experiment_setup
```

**CLI overrides** (take precedence over `calculator_params.toml`):

```bash
# Use a different reagents file
--reagents-file reagents_higher.csv

# Inject fixed-concentration reagents into the conditions CSV at runtime
--fixed-reagents fixed_rxn_concs.csv

# Also overwrite existing columns that match fixed reagents (default: fill missing only)
--fixed-reagents-overwrite-existing
```

---

## Input Files Reference

### `calculator_params.toml`

Must be present in the experiment directory. Auto-detected as `calculator_params.toml`,
`experiment_params.toml`, or `params.toml` (in that priority order).

```toml
# Required
final_rxn_vol_ul = 20.0

# Which conditions CSV to use (default: "samples_final_concs.csv")
samples_final_concs_file = "bounds_with_standards.csv"

# Which reagents file to use (default: "reagents.csv")
reagents_file = "reagents_higher.csv"

# Pipetting settings
pipetting_scalar = 1.1        # overfill fraction (1.1 = 10% extra)
min_pipetting_vol_ul = 1.0    # minimum pipettable volume
conc_decimals = 3             # rounding for concentration columns
vol_decimals = 3              # rounding for volume columns

# Buffer (makes up remaining volume in each well)
buffer_reagent = "Tris HCl"   # must match a row in reagents.csv

# Master mix settings
use_base_master_mix = true
exclude_from_base_reagents = ["vioC", "vioD"]   # titrated individually even if constant

# Replicate expansion (adds a `replicate` column, default 1 = no replication)
# replicates = 3

# Fixed reagents injected into the conditions CSV before calculations
# (equivalent to --fixed-reagents on the CLI)
# fixed_reagents = "fixed_rxn_concs.csv"
# fixed_reagents_overwrite_existing = false
```

### `reagents.csv`

One row per reagent with its stock concentration. Must include the buffer reagent.

```csv
reagent,stock_conc,units
hepes,1500.0,mM
magnesium_acetate,1000.0,mM
...
water,,x
```

For standard curves, also add the pure compound stocks:

```csv
violacein,125.0,uM
deoxyviolacein,125.0,uM
```

### Conditions CSV (output of `doe.py`)

One row per condition. Column names must be `<reagent> <unit>` matching `reagents.csv`.
When mixing DOE and standards, also include a `Type` column (`"doe"` or `"standard"`).

```csv
conditions,Type,vioC uM,vioD uM,...,violacein uM,deoxyviolacein uM
0,doe,0.1,0.8,...,0.0,0.0
...
90,standard,0.0,0.0,...,6.5,0.0
```

---

## Output Files

| File | Description |
|------|-------------|
| `base_master_mix.csv` | Volumes to prepare the pooled base master mix |
| `samples_titration.csv` | Per-well pipetting volumes for each condition |
| `samples_titration_labcraft.csv` | Well-indexed for liquid handlers |
| `reagents.csv` (updated) | `master_mix` row added/updated in place |

### Reading the summary printout

```
base_master_mix_fold=3.4           ← master mix is concentrated 3.4× before dilution
base_master_mix_vol_to_add_ul=17.4 ← µL of master mix per well
base_buffer_vol_per_rxn_ul=11.5   ← buffer in the master mix per well equivalent
max_titration_vol_per_rxn_ul=2.6  ← largest titration volume across all wells
min_titration_vol_per_rxn_ul=0.04 ← smallest titration volume (watch for pipetting limits)
low_component_volumes_below_min_pipetting=0  ← should be 0; non-zero means a stock is too dilute
```

---

## Full Example (LHC)

```bash
# 1. Generate 20 LHC conditions from bounds
PYTHONPATH=src python src/cdk/calculators/doe.py lhc \
  examples/discovery_plate_pipeline/bounds.csv \
  20 \
  --fixed_reagents examples/discovery_plate_pipeline/fixed_final_rxn_concs.csv

# 2. Point calculator_params.toml at the output, then run the pipeline
PYTHONPATH=src python src/cdk/calculators/discovery_plate_concentration_pipeline.py \
  --experiment-dir examples/discovery_plate_pipeline
```

## Full Example (Ratio Sweep)

```bash
# 1. Write ratio_sweep_params.toml in your experiment directory (see above)

# 2. Generate conditions
PYTHONPATH=src python src/cdk/calculators/doe.py ratio-sweep \
  path/to/experiment_setup/ratio_sweep_params.toml

# 3. Set samples_final_concs_file in calculator_params.toml, then run
PYTHONPATH=src python src/cdk/calculators/discovery_plate_concentration_pipeline.py \
  --experiment-dir path/to/experiment_setup
```

## Full Example (Ratio Sweep + Standards)

```bash
# 1. Generate ratio sweep conditions
PYTHONPATH=src python src/cdk/calculators/doe.py ratio-sweep \
  path/to/experiment_setup/ratio_sweep_params.toml

# 2. Generate standards dilution curves
PYTHONPATH=src python src/cdk/calculators/doe.py standards-curve \
  path/to/experiment_setup/standards_params.toml

# 3. Combine into one conditions CSV
PYTHONPATH=src python src/cdk/calculators/doe.py concat \
  path/to/experiment_setup/bounds_<name>.csv \
  path/to/experiment_setup/standards_n12.csv \
  --output path/to/experiment_setup/bounds_with_standards.csv

# 4. Set samples_final_concs_file = "bounds_with_standards.csv" in calculator_params.toml

# 5. Run the pipeline
PYTHONPATH=src python src/cdk/calculators/discovery_plate_concentration_pipeline.py \
  --experiment-dir path/to/experiment_setup
```
