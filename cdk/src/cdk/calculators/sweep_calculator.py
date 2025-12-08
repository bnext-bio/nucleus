from collections import namedtuple
import pandas as pd
import itertools
import os
import csv
import sys
import numpy as np

# Reagent = namedtuple('Reagent', ['plate', 'well', 'stock_conc', 'units', 'description', 'artifact_num'])
#
# REAGENT_INFO = {
#     "hepes": Reagent("reagent_plate", "A1", 1000, "mM", "HEPES", "AR-148"),
#     "potassium_glutamate": Reagent("reagent_plate", "A2", 2500, "mM", "Potassium glutamate", "AR-145"),
#     "magnesium_acetate": Reagent("reagent_plate", "A3", 1000, "mM", "Magnesium acetate", "AR-146"),
#     "ntp": Reagent("reagent_plate", "A4", 100, "mM", "NTP", "AR-728"),
#     "creatine_phosphate": Reagent("reagent_plate", "A5", 1000, "mM", "Creatine phosphate", "AR-702"),
#     "tcep": Reagent("reagent_plate", "A6", 500, "mM", "TCEP", "-"),
#     "folinic_acid": Reagent("reagent_plate", "A7", 5, "mM", "Folinic acid", "AR-161"),
#     "spermidine": Reagent("reagent_plate", "A8", 200, "mM", "Spermidine", "AR-159"),
#     "amino_acid_solution": Reagent("reagent_plate", "A9", 3.25, "mM", "Amino Acid solution", "AR-684"),
#     "trna": Reagent("reagent_plate", "A10", 35, "ug/ul", "tRNA", "AR-730"),
#     "water": Reagent("reagent_plate", "A11", None, "-", "Water", None),
# }
#
# smix_reagents_to_mix = [
#     "hepes", "potassium_glutamate", "magnesium_acetate", "ntp",
#     "creatine_phosphate", "tcep", "folinic_acid", "spermidine",
#     "amino_acid_solution", "water"
# ]
# smix_final_volume = 200  # ul
# smix_fold = 3.3
# final_pure_concs = {
#     'hepes': 50,
#     'potassium_glutamate': 100,
#     'magnesium_acetate': [7.5, 10, 12.5],  # Example sweep
#     'ntp': [1.5, 2, 2.5],                  # Example sweep
#     'creatine_phosphate': 20,
#     'tcep': 1,
#     'folinic_acid': 0.02,
#     'spermidine': 2,
#     'amino_acid_solution': 0.3,
#     'trna': 3.5,
#     'water': None,
#     'neb_sol_a': 1,
#     'neb_sol_b': 1, 'dna': 5, 'rnas_inh': 1,
# }
# reaction_volume = 10  # ul
def check_valid_input_params(REAGENT_INFO, final_pure_concs):
    # ============================================================================
    # ASSERTION 1: All stock concentrations must be positive
    # ============================================================================
    print("\n[1] Checking stock concentrations are positive...")
    for reagent_name, reagent_info in REAGENT_INFO.items():
        if reagent_info.stock_conc is not None:
            assert reagent_info.stock_conc > 0, \
                f"ERROR: {reagent_name} has non-positive stock concentration: {reagent_info.stock_conc}"
    print("✓ All stock concentrations are positive")

    # ============================================================================
    # ASSERTION 2: All target concentrations must be positive (if specified)
    # ============================================================================
    print("\n[2] Checking target concentrations are positive...")
    for reagent, conc in final_pure_concs.items():
        if conc is not None:
            if isinstance(conc, list):
                for c in conc:
                    assert c > 0, f"ERROR: {reagent} has non-positive target concentration: {c}"
            else:
                assert conc > 0, f"ERROR: {reagent} has non-positive target concentration: {conc}"
    print("✓ All target concentrations are positive")

def get_sweep_reagents(final_pure_concs):
    # Identify sweep reagents
    sweep_reagents = []
    sweep_values = []
    for reagent, conc in final_pure_concs.items():
        if isinstance(conc, list) and len(conc) > 1:
            sweep_reagents.append(reagent)
            sweep_values.append(conc)
    return sweep_reagents, sweep_values

def calc_base_conc_for_sweep(smix_reagents_to_mix, final_pure_concs):
    # Determine minimal concentration for each sweep reagent (for smix)
    smix_base_concs = {}
    for r in smix_reagents_to_mix:
        if r in final_pure_concs:
            conc = final_pure_concs[r]
            if isinstance(conc, list):
                smix_base_concs[r] = min(conc)  # Use minimal concentration for smix
            else:
                smix_base_concs[r] = conc
        else:
            smix_base_concs[r] = None
    return smix_base_concs

def calc_smix_conc(smix_base_concs, smix_reagents_to_mix, smix_fold, reagents, smix_final_volume, reaction_volume):
    # Calculate smix composition at 3.3x concentration (using minimal values)
    smix_target_conc = {}
    smix_vols_needed = {}

    for r in smix_reagents_to_mix:
        pure_conc = smix_base_concs[r]
        if pure_conc is not None:
            smix_target_conc[r] = pure_conc * smix_fold
            stock = reagents[r].stock_conc
            smix_vols_needed[r] = round((smix_target_conc[r] * smix_final_volume) / stock, 3)
        else:
            smix_target_conc[r] = None
            smix_vols_needed[r] = None

    # ============================================================================
    # ASSERTION 4: All smix volumes must be positive
    # ============================================================================
    print("\n[4] Checking all smix volumes are positive...")
    for reagent, vol in smix_vols_needed.items():
        if vol is not None:
            assert vol > 0, f"ERROR: {reagent} has non-positive smix volume: {vol} ul"
    print("✓ All smix volumes are positive")

    # Calculate water volume
    total_non_water_vol = sum(v for r, v in smix_vols_needed.items() if r != 'water' and v is not None)
    water_vol = round(smix_final_volume - total_non_water_vol, 3)
    smix_vols_needed['water'] = water_vol

    # ============================================================================
    # ASSERTION 5: Smix total volume must equal target volume
    # ============================================================================
    print("\n[5] Checking smix total volume equals target...")
    smix_total = sum(v for v in smix_vols_needed.values() if v is not None)
    assert abs(smix_total - smix_final_volume) < 0.01, \
        f"ERROR: Smix total volume ({smix_total} ul) != target volume ({smix_final_volume} ul)"
    print(f"✓ Smix total volume = {smix_total:.3f} ul (target: {smix_final_volume} ul)")

    # ============================================================================
    # ASSERTION 6: Water volume must be positive
    # ============================================================================
    print("\n[6] Checking water volume is positive...")
    assert water_vol > 0, f"ERROR: Water volume is non-positive: {water_vol} ul (reagents overfilled!)"
    print(f"✓ Water volume = {water_vol:.3f} ul")

    # Calculate smix volume needed per reaction
    smix_vol_per_reaction = round(reaction_volume / smix_fold, 3)  # ul
    return smix_target_conc, smix_vols_needed, smix_vol_per_reaction

def check_reaction_conditions(reaction_rows, sweep_reagents,reaction_volume, smix_vol_per_reaction,
                              final_pure_concs, reagents, warnings=[]):
    # ============================================================================
    # ASSERTION 7: All extra volumes must be non-negative
    # ============================================================================
    # print("\n[7] Checking all extra volumes are non-negative...")
    # for row in reaction_rows:
    #     for reagent in sweep_reagents:
    #         extra_vol = row[f'{reagent}_extra_ul']
    #         assert extra_vol >= 0, \
    #             f"ERROR: {row['condition']} has negative extra volume for {reagent}: {extra_vol} ul"
    # print("✓ All extra volumes are non-negative")

    # ============================================================================
    # ASSERTION 8: Total volume (smix + extras) must not exceed reaction volume
    # ============================================================================
    # print("\n[8] Checking total volumes don't exceed reaction volume...")
    # max_extra_vol_allowed = reaction_volume - smix_vol_per_reaction
    # for row in reaction_rows:
    #     total_check = row['total_vol_check']
    #     if total_check > reaction_volume:
    #         warnings.append(
    #             f"WARNING: {row['condition']} total volume ({total_check} ul) exceeds reaction volume ({reaction_volume} ul)")
    #     assert row['total_extra_vol'] <= max_extra_vol_allowed, \
    #         f"ERROR: {row['condition']} extra volumes ({row['total_extra_vol']} ul) exceed available space ({max_extra_vol_allowed} ul)"
    # print(f"✓ All extra volumes fit within reaction volume (max allowed: {max_extra_vol_allowed:.3f} ul)")

    # ============================================================================
    # ASSERTION 9: Smix volume must be less than reaction volume
    # ============================================================================
    print("\n[9] Checking smix volume leaves room for titrations...")
    assert smix_vol_per_reaction < reaction_volume, \
        f"ERROR: Smix volume ({smix_vol_per_reaction} ul) >= reaction volume ({reaction_volume} ul)"
    print(f"✓ Smix volume ({smix_vol_per_reaction} ul) < reaction volume ({reaction_volume} ul)")

    # ============================================================================
    # ASSERTION 10: Target concentrations achievable with stock concentrations
    # ============================================================================
    print("\n[10] Checking target concentrations are achievable...")
    for reagent, conc in final_pure_concs.items():
        if conc is not None and reagent in reagents:
            stock_conc = reagents[reagent].stock_conc
            if stock_conc is not None:
                if isinstance(conc, list):
                    max_target = max(conc)
                else:
                    max_target = conc
                max_achievable = stock_conc
                assert max_target <= max_achievable, \
                    f"ERROR: {reagent} target ({max_target} mM) exceeds stock concentration ({stock_conc} mM)"
    print("✓ All target concentrations are achievable with available stocks")

    # Display warnings if any
    if warnings:
        print("\n" + "!" * 60)
        print("WARNINGS:")
        for warning in warnings:
            print(f"  {warning}")
        print("!" * 60)

    # return warnings
def add_fixed_volumes_to_df(reaction_rows, sweep_reagents, smix_reagents_to_mix, final_pure_concs, reagents,
                            reaction_volume):
    # Build main reaction conditions DataFrame
    reaction_conditions = pd.DataFrame(reaction_rows)

    non_smix_reagents = [
        r for r in final_pure_concs.keys()
        if r not in smix_reagents_to_mix and r != 'water'
    ]
    # Add vectorized columns for all fixed reagents (NOT in smix, NOT sweep)
    fixed_reagents = [r for r in non_smix_reagents if
                      not (r in sweep_reagents and r in smix_reagents_to_mix) and r != 'water' and r in reagents]

    for reagent in fixed_reagents:
        conc = final_pure_concs[reagent]
        stock_conc = reagents[reagent].stock_conc
        if stock_conc and conc is not None:
            vol_needed = round((conc * reaction_volume) / stock_conc, 4)
            reaction_conditions[f'{reagent}_vol_ul'] = vol_needed


    # After all reagent columns are added, compute the sum of fixed additions for each row
    reaction_conditions['fixed_total_ul'] = reaction_conditions[[f'{r}_vol_ul' for r in fixed_reagents]].sum(axis=1)

    # Now calculate water column in vectorized fashion
    # For each row: water = reaction_volume - (smix + all titrations + all fixed)
    tit_columns = [col for col in reaction_conditions.columns if col.endswith('_extra_ul')]
    reaction_conditions['titration_sum_ul'] = reaction_conditions[tit_columns].sum(axis=1) if tit_columns else 0
    reaction_conditions['water_vol_ul'] = (
            reaction_volume
            - reaction_conditions['bnext_sms_vol_ul']
            - reaction_conditions['titration_sum_ul']
            - reaction_conditions['fixed_total_ul']
    ).round(4).clip(lower=0.0)

    # Optionally compute final volume check
    reaction_conditions['final_vol_check'] = (
            reaction_conditions['bnext_sms_vol_ul']
            + reaction_conditions['titration_sum_ul']
            + reaction_conditions['fixed_total_ul']
            + reaction_conditions['water_vol_ul']
    ).round(4)
    assert reaction_conditions['final_vol_check'].any() < reaction_volume, 'OT appropriate print statement. Uncomment next line instead'
    # assert reaction_conditions['final_vol_check'].any() < reaction_volume, f'Lower sweep reagent concentrations -- not enough dead volume for reaction condition {reaction_conditions[reaction_conditions['final_vol_check'] > reaction_volume]}'
    return reaction_conditions

def find_suitable_reagent(base_name, reagents, target_conc, reaction_volume, min_vol,
                          used_volume):
    """
    Given a sweep reagent, find the best stock (original or diluted)
    that gives >= min_vol. Creates a new aliquot if none exist.
    """
    # if sweep_reagent is not None:
    #     base_name = sweep_reagent["name"]
    print(f'Finding suitable available concentrations of {base_name} for {target_conc} conc with minimum pipetting volume of {min_vol}ul')
    # Gather all available versions of this reagent (base + aliquots)
    candidate_names = [r for r in reagents if r.startswith(base_name)]
    candidates = [(name, reagents[name]) for name in candidate_names]
    # log(f"{base_name} has {candidates} candidates")
    # Try each candidate stock concentration
    for name, reagent in sorted(candidates, key=lambda x: -x[1].stock_conc):
        vol = (reaction_volume * target_conc) / reagent.stock_conc
        water_to_add = reaction_volume - used_volume - vol
        print(f'    If using {name} will need vol {vol} + water {water_to_add}')
        if (vol >= min_vol) and (water_to_add >= 0):
            return reagent, vol, name  # found a usable stock

    # If none suitable, determine a suitable diluted stock factor
    parent = reagents[base_name] #TODO: choose most appropriate concentration dont just take base
    # Dilute just enough so the required volume is >= min_vol
    child_conc = np.ceil(reaction_volume * target_conc) / min_vol
    dilution_factor = parent.stock_conc / child_conc
    assert dilution_factor > 1, f"No available stock to be made from {parent.stock_conc} {parent.units} for {target_conc} {parent.units}"

    # new_reagent = create_and_register_diluted_reagent(reagents, base_name, dilution_factor, protocol, log, pipette)
    vol = (reaction_volume * target_conc) / child_conc
    assert False, (f'Failed to find a suitable reagent for {target_conc} {parent.units} within minimum pipetting bounds'
                   f'\n Suggest adding a {base_name} reagent w/ dilution factor {dilution_factor} for {child_conc} {parent.units}')
    # return reagent, vol
    # assert np.all(vol > 0), f"Impossible mixture for concs {target_conc} vol is {vol}"
    # print(f'Make a diluted stock of {child_conc} {parent.units} with dilution_factor: {dilution_factor} for {target_conc} {parent.units}')


def generate_sweep_combinations(sweep_reagents, sweep_values, smix_vol_per_reaction, smix_base_concs,
                                reagents, reaction_volume):
    # Generate all combinations of sweep reagents
    if sweep_reagents:
        combinations = list(itertools.product(*sweep_values))
    else:
        combinations = [tuple()]

    # For each combination, calculate extra volumes needed
    reaction_rows = []

    for combo_idx, combo in enumerate(combinations):
        row = {'condition': f'Condition_{combo_idx + 1}', 'bnext_sms_vol_ul': smix_vol_per_reaction}
        # Add target concentrations for sweep reagents
        for i, reagent in enumerate(sweep_reagents):
            row[f'{reagent}_target_mM'] = combo[i]
        # Calculate extra volume needed for each sweep reagent
        for i, reagent in enumerate(sweep_reagents):
            target_conc = combo[i]
            base_conc = smix_base_concs[reagent]  # Concentration provided by smix
            extra_conc_needed = target_conc - base_conc
            if extra_conc_needed > 0:
                best_reagent, extra_vol, name = (
                    find_suitable_reagent(reagent, reagents, extra_conc_needed,
                                          reaction_volume, min_vol=0.3, used_volume=0)) #currently not tracking everything in well yet
                # stock_conc = reagents[reagent].stock_conc
                # extra_vol = round((extra_conc_needed * reaction_volume) / stock_conc, 4)
                row[f'{name}_extra_ul'] = extra_vol
            else:
                row[f'{reagent}_extra_ul'] = 0.0
        reaction_rows.append(row)
        print(f'Sucessfully generated {combo} concs for {sweep_reagents} sweep reagents')
    return reaction_rows

def save_table(table, file_path):
    # def save_tracked_wells(self, file_path, log):
    """Save wells to CSV if local, otherwise log to comments on robot."""
    if os.environ.get("OT_SIMULATE_CSV"):
        print(f'Tables saved to {file_path}')
        table.to_csv(file_path, index=False)
    else:
        print(f'Attempted to save to {file_path},'
            f'But the environment variable "OT_SIMULATE_CSV" is not set.'                
            f'To save the file on this machine, run: export OT_SIMULATE_CSV=1'
            f'Do not set on OT since it is a write only machine')

def create_output_tables(smix_reagents_to_mix, smix_base_concs, smix_target_conc,
                         smix_vols_needed, reaction_conditions, save_tables=True):
    # Create output tables
    smix_table = pd.DataFrame({
        'reagent': smix_reagents_to_mix,
        'smix_base_conc_mM': [smix_base_concs[r] for r in smix_reagents_to_mix],
        'smix_target_conc_mM': [smix_target_conc[r] for r in smix_reagents_to_mix],
        'bnext_sms_vol_ul': [smix_vols_needed[r] for r in smix_reagents_to_mix],
    })

    # reaction_conditions = pd.DataFrame(reaction_rows)
    if save_tables:
        # Save outputs
        # smix_table.to_csv("smix_volumes_sweep.csv", index=False)
        # reaction_conditions.to_csv("reaction_conditions_sweep.csv", index=False)
        main_dir = f'/Volumes/bnext/experiments/ot_smix_params/energy_mix_'
        save_table(smix_table, main_dir + "smix_volumes_sweep.csv")
        save_table(reaction_conditions, main_dir + "reaction_conditions_sweep.csv")
    return smix_table, reaction_conditions

def generate_smix_and_sweep_tables(smix_reagents_to_mix, smix_fold, smix_final_volume,
                                   final_pure_concs, reagents, reaction_volume):
    smix_base_concs= calc_base_conc_for_sweep(smix_reagents_to_mix, final_pure_concs)
    smix_target_conc, smix_vols_needed, smix_vol_per_reaction = calc_smix_conc(smix_base_concs, smix_reagents_to_mix,
                                                                               smix_fold, reagents, smix_final_volume,
                                                                               reaction_volume)

    sweep_reagents, sweep_values = get_sweep_reagents(final_pure_concs)
    reaction_rows = generate_sweep_combinations(sweep_reagents, sweep_values, smix_vol_per_reaction, smix_base_concs,
                                                reagents, reaction_volume
                                                )
    reaction_conditions = add_fixed_volumes_to_df(reaction_rows, sweep_reagents, smix_reagents_to_mix, final_pure_concs, reagents,
                            reaction_volume)
    smix_table, reaction_conditions = create_output_tables(smix_reagents_to_mix, smix_base_concs, smix_target_conc,
                                                        smix_vols_needed, reaction_conditions, save_tables=True)
    check_reaction_conditions(reaction_rows, sweep_reagents, reaction_volume, smix_vol_per_reaction,
                              final_pure_concs, reagents, warnings=[])
    return smix_table, reaction_conditions

# smix_table, reaction_conditions  = generate_smix_and_sweep_tables(smix_reagents_to_mix, smix_fold, smix_final_volume, final_pure_concs, reagents, reaction_volume)