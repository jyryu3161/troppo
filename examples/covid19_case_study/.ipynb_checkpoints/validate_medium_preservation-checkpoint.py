#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Medium Preservation Validation
==============================
Validates whether tINIT reconstruction preserved the original model's
exchange reactions, uptake bounds, and medium conditions.

Outputs: covid19_tinit_results/medium_validation_report.txt
         covid19_tinit_results/medium_validation_details.csv
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(BASE_DIR, "HumanGEM_Consistent_COVID19_HAM.xml")
MODELS_DIR = os.path.join(BASE_DIR, "covid19_tinit_results", "models")
OUTPUT_DIR = os.path.join(BASE_DIR, "covid19_tinit_results")
OMICS_PATH = os.path.join(BASE_DIR, "Desai-GTEx_ensembl.csv")

sys.path.append(os.path.join(BASE_DIR, "../../src"))

from cobra.io import read_sbml_model

# ══════════════════════════════════════════════════════════════════
# 1. Load original model
# ══════════════════════════════════════════════════════════════════
print("=" * 80)
print("Medium Preservation Validation Report")
print("=" * 80)

print("\n[1] Loading original model...")
original = read_sbml_model(MODEL_PATH)
print(f"  Original: {len(original.reactions)} rxns, {len(original.metabolites)} mets, {len(original.genes)} genes")

# Identify all exchange reactions in original
orig_exchanges = {}
for rxn in original.reactions:
    if rxn.boundary:
        orig_exchanges[rxn.id] = {
            'name': rxn.name,
            'lb': rxn.lower_bound,
            'ub': rxn.upper_bound,
            'metabolites': ', '.join(m.id for m in rxn.metabolites),
            'active_uptake': rxn.lower_bound < 0,
        }

print(f"  Total exchange reactions: {len(orig_exchanges)}")
active_uptake = sum(1 for v in orig_exchanges.values() if v['active_uptake'])
print(f"  Active uptake (LB < 0): {active_uptake}")

# ══════════════════════════════════════════════════════════════════
# 2. Define essential nutrient categories
# ══════════════════════════════════════════════════════════════════
# Key nutrients for Ham's F12 + DMEM medium and human metabolism
essential_nutrients = {
    # Essential amino acids
    'HMR_9036': 'L-Histidine',
    'HMR_9038': 'L-Isoleucine',
    'HMR_9039': 'L-Leucine',
    'HMR_9040': 'L-Lysine',
    'HMR_9041': 'L-Methionine',
    'HMR_9043': 'L-Phenylalanine',
    'HMR_9044': 'L-Threonine',
    'HMR_9045': 'L-Tryptophan',
    'HMR_9046': 'L-Valine',
    # Non-essential amino acids (medium components)
    'HMR_9033': 'L-Alanine',
    'HMR_9035': 'L-Arginine',
    'HMR_9032': 'L-Asparagine',
    'HMR_9037': 'L-Glutamine',
    'HMR_9042': 'Glycine',
    'HMR_9047': 'L-Proline',
    'HMR_9063': 'L-Serine',
    'HMR_9064': 'L-Tyrosine',
    'HMR_9065': 'L-Cysteine',
    # Core carbon/energy
    'HMR_9034': 'D-Glucose',
    'HMR_9048': 'O2 (oxygen)',
    'HMR_9058': 'CO2 (carbon dioxide)',
    # Vitamins
    'HMR_9109': 'Thiamine (B1)',
    'HMR_9110': 'Riboflavin (B2)',
    'HMR_9143': 'Niacinamide (B3/NAM)',
    'HMR_9144': 'Pantothenate (B5)',
    'HMR_9145': 'Pyridoxine (B6)',
    'HMR_9146': 'Folate (B9)',
    'HMR_9151': 'Choline',
    'HMR_9153': 'Inositol',
    # Minerals & ions
    'HMR_9078': 'Phosphate (Pi)',
    'HMR_9079': 'Sulfate (SO4)',
    'HMR_9072': 'Iron (Fe2+)',
    'HMR_9073': 'Iron (Fe3+)',
    'HMR_9074': 'Potassium (K+)',
    'HMR_9075': 'Sodium (Na+)',
    'HMR_9076': 'Calcium (Ca2+)',
    # Water
    'HMR_9135': 'H2O (water)',
}

# ══════════════════════════════════════════════════════════════════
# 3. Load reconstructed models & compare
# ══════════════════════════════════════════════════════════════════
print("\n[2] Loading reconstructed models and comparing...")

sample_names = [
    'Lung_Healthy', 'Lung_COVID19',
    'Heart_Healthy', 'Heart_COVID19',
    'Liver_Healthy', 'Liver_COVID19',
    'Kidney_Healthy', 'Kidney_COVID19',
    'Intestine_Healthy', 'Intestine_COVID19',
]

all_comparison_rows = []
missing_by_sample = {}
bounds_changed_by_sample = {}
biomass_results = {}

for sample in sample_names:
    xml_path = os.path.join(MODELS_DIR, f"{sample}_model.xml")
    if not os.path.exists(xml_path):
        print(f"  {sample}: model file not found, skipping")
        continue

    recon = read_sbml_model(xml_path)
    recon_rxn_ids = set(r.id for r in recon.reactions)
    recon_exchange_map = {}
    for rxn in recon.reactions:
        if rxn.boundary:
            recon_exchange_map[rxn.id] = {'lb': rxn.lower_bound, 'ub': rxn.upper_bound}

    missing = []
    bounds_changed = []

    for rxn_id, info in orig_exchanges.items():
        row = {
            'sample': sample,
            'rxn_id': rxn_id,
            'rxn_name': info['name'],
            'metabolites': info['metabolites'],
            'orig_lb': info['lb'],
            'orig_ub': info['ub'],
            'orig_active_uptake': info['active_uptake'],
            'is_essential_nutrient': rxn_id in essential_nutrients,
            'nutrient_label': essential_nutrients.get(rxn_id, ''),
        }

        if rxn_id not in recon_rxn_ids:
            row['status'] = 'REMOVED'
            row['recon_lb'] = None
            row['recon_ub'] = None
            missing.append(rxn_id)
        elif rxn_id in recon_exchange_map:
            r_info = recon_exchange_map[rxn_id]
            row['recon_lb'] = r_info['lb']
            row['recon_ub'] = r_info['ub']
            if abs(info['lb'] - r_info['lb']) > 1e-6 or abs(info['ub'] - r_info['ub']) > 1e-6:
                row['status'] = 'BOUNDS_CHANGED'
                bounds_changed.append(rxn_id)
            else:
                row['status'] = 'PRESERVED'
        else:
            # Reaction exists but not as exchange (unlikely)
            row['status'] = 'PRESERVED_NON_BOUNDARY'
            row['recon_lb'] = None
            row['recon_ub'] = None

        all_comparison_rows.append(row)

    missing_by_sample[sample] = missing
    bounds_changed_by_sample[sample] = bounds_changed

    # ── Biomass FBA ──
    recon.objective = 'biomass_human'
    sol = recon.optimize()
    bm = sol.objective_value if sol.objective_value is not None else 0.0
    biomass_results[sample] = bm

    n_recon_ex = len(recon_exchange_map)
    print(f"  {sample}: {len(recon.reactions)} rxns, "
          f"exchanges={n_recon_ex}, "
          f"removed={len(missing)}, "
          f"bounds_changed={len(bounds_changed)}, "
          f"biomass={bm:.6f}")

# ══════════════════════════════════════════════════════════════════
# 4. Build detailed comparison DataFrame
# ══════════════════════════════════════════════════════════════════
comp_df = pd.DataFrame(all_comparison_rows)
detail_csv = os.path.join(OUTPUT_DIR, "medium_validation_details.csv")
comp_df.to_csv(detail_csv, index=False)
print(f"\n  Detailed comparison saved: {detail_csv}")

# ══════════════════════════════════════════════════════════════════
# 5. Functional Uptake Test (FBA)
# ══════════════════════════════════════════════════════════════════
print("\n[3] Functional Uptake Test (FBA)...")

uptake_test_results = []
test_nutrients = {k: v for k, v in essential_nutrients.items() if k in orig_exchanges}

# Test on original model first
print("  Testing original model...")
original.objective = 'biomass_human'
orig_sol = original.optimize()
orig_bm = orig_sol.objective_value if orig_sol.objective_value is not None else 0.0
print(f"  Original biomass: {orig_bm:.6f}")

for rxn_id, label in test_nutrients.items():
    if rxn_id not in set(r.id for r in original.reactions):
        continue
    rxn = original.reactions.get_by_id(rxn_id)
    orig_lb = rxn.lower_bound
    # Test: block this nutrient
    with original:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
        sol = original.optimize()
        blocked_bm = sol.objective_value if sol.objective_value is not None else 0.0
    uptake_test_results.append({
        'sample': 'ORIGINAL',
        'rxn_id': rxn_id,
        'nutrient': label,
        'orig_lb': orig_lb,
        'biomass_normal': orig_bm,
        'biomass_blocked': blocked_bm,
        'essential_for_growth': blocked_bm < 0.01 * orig_bm if orig_bm > 1e-6 else None,
    })

# Test on each reconstructed model
for sample in sample_names:
    xml_path = os.path.join(MODELS_DIR, f"{sample}_model.xml")
    if not os.path.exists(xml_path):
        continue

    recon = read_sbml_model(xml_path)
    recon.objective = 'biomass_human'
    recon_rxn_ids = set(r.id for r in recon.reactions)

    sol = recon.optimize()
    normal_bm = sol.objective_value if sol.objective_value is not None else 0.0

    for rxn_id, label in test_nutrients.items():
        if rxn_id not in recon_rxn_ids:
            uptake_test_results.append({
                'sample': sample,
                'rxn_id': rxn_id,
                'nutrient': label,
                'orig_lb': orig_exchanges.get(rxn_id, {}).get('lb', 0),
                'biomass_normal': normal_bm,
                'biomass_blocked': None,
                'essential_for_growth': None,
                'note': 'REACTION_ABSENT',
            })
        else:
            rxn = recon.reactions.get_by_id(rxn_id)
            with recon:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
                sol2 = recon.optimize()
                blocked_bm = sol2.objective_value if sol2.objective_value is not None else 0.0
            uptake_test_results.append({
                'sample': sample,
                'rxn_id': rxn_id,
                'nutrient': label,
                'orig_lb': orig_exchanges.get(rxn_id, {}).get('lb', 0),
                'biomass_normal': normal_bm,
                'biomass_blocked': blocked_bm,
                'essential_for_growth': blocked_bm < 0.01 * normal_bm if normal_bm > 1e-6 else None,
            })

    print(f"  {sample}: uptake tests done (biomass={normal_bm:.6f})")

uptake_df = pd.DataFrame(uptake_test_results)
uptake_csv = os.path.join(OUTPUT_DIR, "uptake_test_results.csv")
uptake_df.to_csv(uptake_csv, index=False)
print(f"  Uptake test results saved: {uptake_csv}")

# ══════════════════════════════════════════════════════════════════
# 6. Gene Expression vs Removed Exchange Reactions
# ══════════════════════════════════════════════════════════════════
print("\n[4] Gene expression analysis for removed exchange reactions...")

omics_data = pd.read_csv(OMICS_PATH, index_col=0)

# Load thresholded data
from troppo.omics import GeneLevelThresholding
glt = GeneLevelThresholding(
    omics_dataframe=omics_data,
    thresholding_strat='local t2',
    global_threshold_lower=2,
    global_threshold_upper=3,
    local_threshold=2,
)
filtered_df = glt.apply_thresholding_filter()

conflict_analysis = []
filter_suffix = '_local_t2_2_3_2'

for sample in sample_names:
    missing = missing_by_sample.get(sample, [])
    filtered_key = sample + filter_suffix

    for rxn_id in missing:
        if rxn_id not in set(r.id for r in original.reactions):
            continue
        rxn = original.reactions.get_by_id(rxn_id)
        gpr = rxn.gene_reaction_rule

        row = {
            'sample': sample,
            'rxn_id': rxn_id,
            'rxn_name': rxn.name,
            'gpr': gpr if gpr else 'NO_GPR',
            'orig_lb': rxn.lower_bound,
            'orig_ub': rxn.upper_bound,
            'was_active_uptake': rxn.lower_bound < 0,
            'is_essential_nutrient': rxn_id in essential_nutrients,
            'nutrient_label': essential_nutrients.get(rxn_id, ''),
        }

        # Check gene expression
        if gpr:
            import re
            gene_ids = re.findall(r'ENSG\d+', gpr)
            gene_scores = {}
            for g in gene_ids:
                if g in filtered_df.columns and filtered_key in filtered_df.index:
                    gene_scores[g] = filtered_df.loc[filtered_key, g]
                elif g in omics_data.columns and sample in omics_data.index:
                    gene_scores[g] = omics_data.loc[sample, g]
            row['genes'] = '; '.join(gene_ids)
            row['gene_scores'] = '; '.join(f"{g}={v:.2f}" for g, v in gene_scores.items())
            row['min_score'] = min(gene_scores.values()) if gene_scores else None
            row['max_score'] = max(gene_scores.values()) if gene_scores else None
            row['all_below_threshold'] = all(v < 0 for v in gene_scores.values()) if gene_scores else None
            row['removal_cause'] = 'LOW_EXPRESSION' if row.get('all_below_threshold') else 'MIXED/OTHER'
        else:
            row['genes'] = ''
            row['gene_scores'] = ''
            row['min_score'] = None
            row['max_score'] = None
            row['all_below_threshold'] = None
            row['removal_cause'] = 'NO_GPR (score=0 in tINIT)'

        conflict_analysis.append(row)

conflict_df = pd.DataFrame(conflict_analysis)
conflict_csv = os.path.join(OUTPUT_DIR, "expression_conflict_analysis.csv")
if not conflict_df.empty:
    conflict_df.to_csv(conflict_csv, index=False)
    print(f"  Expression conflict analysis saved: {conflict_csv}")
else:
    print("  No removed exchange reactions to analyze.")

# ══════════════════════════════════════════════════════════════════
# 7. Generate Summary Report
# ══════════════════════════════════════════════════════════════════
print("\n[5] Generating summary report...")

report_path = os.path.join(OUTPUT_DIR, "medium_validation_report.txt")
with open(report_path, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("Medium Preservation Validation Report\n")
    f.write("=" * 80 + "\n\n")

    # ── Section 1: Overview ──
    f.write("1. OVERVIEW\n")
    f.write("-" * 40 + "\n")
    f.write(f"Original model: {original.id}\n")
    f.write(f"  Reactions: {len(original.reactions)}\n")
    f.write(f"  Exchange reactions: {len(orig_exchanges)}\n")
    f.write(f"  Active uptake (LB<0): {active_uptake}\n")
    f.write(f"  Original biomass flux: {orig_bm:.6f}\n\n")

    # ── Section 2: Exchange reaction preservation ──
    f.write("\n2. EXCHANGE REACTION PRESERVATION (per sample)\n")
    f.write("-" * 40 + "\n")
    f.write(f"{'Sample':<25s} {'Total_Ex':>8s} {'Preserved':>10s} {'Removed':>8s} {'Bounds_Chg':>11s} {'Biomass':>10s}\n")
    for sample in sample_names:
        sample_rows = comp_df[comp_df['sample'] == sample]
        n_preserved = (sample_rows['status'] == 'PRESERVED').sum()
        n_removed = (sample_rows['status'] == 'REMOVED').sum()
        n_changed = (sample_rows['status'] == 'BOUNDS_CHANGED').sum()
        bm = biomass_results.get(sample, 0.0)
        f.write(f"{sample:<25s} {len(orig_exchanges):>8d} {n_preserved:>10d} {n_removed:>8d} {n_changed:>11d} {bm:>10.6f}\n")

    # ── Section 3: Removed exchange reactions (essential nutrients) ──
    f.write("\n\n3. REMOVED ESSENTIAL NUTRIENT EXCHANGES\n")
    f.write("-" * 40 + "\n")
    essential_removed = comp_df[(comp_df['status'] == 'REMOVED') & (comp_df['is_essential_nutrient'])]
    if essential_removed.empty:
        f.write("  None - all essential nutrient exchanges preserved.\n")
    else:
        f.write(f"  Total essential nutrient exchanges removed: {len(essential_removed)}\n\n")
        for _, row in essential_removed.iterrows():
            f.write(f"  {row['sample']:<25s} {row['rxn_id']:<15s} {row['nutrient_label']:<25s} "
                    f"LB={row['orig_lb']:>8.1f}\n")

    # ── Section 4: All removed exchange reactions (full list) ──
    f.write("\n\n4. ALL REMOVED EXCHANGE REACTIONS (MAR ID list)\n")
    f.write("-" * 40 + "\n")
    removed_all = comp_df[comp_df['status'] == 'REMOVED']
    if removed_all.empty:
        f.write("  None - all exchange reactions preserved.\n")
    else:
        # Group by rxn_id to show which samples lost each reaction
        removed_pivot = removed_all.groupby('rxn_id').agg({
            'rxn_name': 'first',
            'orig_lb': 'first',
            'orig_ub': 'first',
            'orig_active_uptake': 'first',
            'is_essential_nutrient': 'first',
            'sample': lambda x: ', '.join(sorted(x)),
        }).reset_index()
        removed_pivot = removed_pivot.sort_values('orig_active_uptake', ascending=False)

        f.write(f"  Unique removed exchange reactions: {len(removed_pivot)}\n\n")
        f.write(f"  {'MAR_ID':<15s} {'Name':<35s} {'LB':>8s} {'UB':>8s} {'Uptake':>7s} {'Essential':>9s} Removed_In\n")
        f.write("  " + "-" * 120 + "\n")
        for _, row in removed_pivot.iterrows():
            f.write(f"  {row['rxn_id']:<15s} {str(row['rxn_name'])[:34]:<35s} "
                    f"{row['orig_lb']:>8.1f} {row['orig_ub']:>8.1f} "
                    f"{'YES' if row['orig_active_uptake'] else 'no':>7s} "
                    f"{'YES' if row['is_essential_nutrient'] else '':>9s} "
                    f"{row['sample']}\n")

    # ── Section 5: Bounds changed ──
    f.write("\n\n5. BOUNDS CHANGED (before vs after)\n")
    f.write("-" * 40 + "\n")
    changed = comp_df[comp_df['status'] == 'BOUNDS_CHANGED']
    if changed.empty:
        f.write("  None - all preserved reactions retain original bounds.\n")
    else:
        f.write(f"  Total reactions with changed bounds: {len(changed)}\n\n")
        f.write(f"  {'Sample':<25s} {'MAR_ID':<15s} {'Orig_LB':>8s} {'Recon_LB':>9s} {'Orig_UB':>8s} {'Recon_UB':>9s}\n")
        f.write("  " + "-" * 80 + "\n")
        for _, row in changed.iterrows():
            f.write(f"  {row['sample']:<25s} {row['rxn_id']:<15s} "
                    f"{row['orig_lb']:>8.1f} {row['recon_lb']:>9.1f} "
                    f"{row['orig_ub']:>8.1f} {row['recon_ub']:>9.1f}\n")

    # ── Section 6: Biomass impact ──
    f.write("\n\n6. BIOMASS IMPACT ASSESSMENT (MAR13082 / biomass_human)\n")
    f.write("-" * 40 + "\n")
    f.write(f"  Original model biomass: {orig_bm:.6f}\n\n")
    f.write(f"  {'Sample':<25s} {'Biomass':>10s} {'vs_Original':>12s} {'Status':>10s}\n")
    f.write("  " + "-" * 60 + "\n")
    for sample in sample_names:
        bm = biomass_results.get(sample, 0.0)
        pct = (bm / orig_bm * 100) if orig_bm > 1e-6 else 0.0
        status = 'OK' if bm > 1e-6 else 'ZERO'
        f.write(f"  {sample:<25s} {bm:>10.6f} {pct:>11.1f}% {status:>10s}\n")

    # ── Section 7: Uptake essentiality ──
    f.write("\n\n7. FUNCTIONAL UPTAKE TEST (essential nutrients)\n")
    f.write("-" * 40 + "\n")
    f.write("  Blocking each nutrient individually and measuring biomass impact.\n\n")

    # Original model results
    orig_uptake = uptake_df[uptake_df['sample'] == 'ORIGINAL']
    f.write("  [Original Model]\n")
    f.write(f"  {'Nutrient':<30s} {'LB':>6s} {'Normal_BM':>10s} {'Blocked_BM':>11s} {'Essential':>10s}\n")
    f.write("  " + "-" * 70 + "\n")
    for _, row in orig_uptake.iterrows():
        ess = 'YES' if row.get('essential_for_growth') else 'no'
        blocked = f"{row['biomass_blocked']:.4f}" if row['biomass_blocked'] is not None else 'N/A'
        f.write(f"  {row['nutrient']:<30s} {row['orig_lb']:>6.0f} {row['biomass_normal']:>10.4f} "
                f"{blocked:>11s} {ess:>10s}\n")

    # Per-sample results (summary)
    f.write("\n  [Per-Sample Summary: nutrients absent from reconstructed model]\n")
    for sample in sample_names:
        sample_uptake = uptake_df[(uptake_df['sample'] == sample)]
        absent = sample_uptake[sample_uptake.get('note', pd.Series()) == 'REACTION_ABSENT'] if 'note' in sample_uptake.columns else pd.DataFrame()
        if not absent.empty:
            f.write(f"\n  {sample}:\n")
            for _, row in absent.iterrows():
                f.write(f"    {row['rxn_id']:<15s} {row['nutrient']:<30s} (LB={row['orig_lb']:.0f}) ABSENT\n")

    # ── Section 8: Expression conflict analysis ──
    f.write("\n\n8. EXPRESSION vs MEDIUM CONFLICT ANALYSIS\n")
    f.write("-" * 40 + "\n")
    f.write("  Removed uptake reactions and their gene expression scores.\n")
    f.write("  'LOW_EXPRESSION': all associated genes below threshold (score < 0)\n")
    f.write("  'NO_GPR': reaction has no gene association (tINIT score = 0)\n\n")

    if not conflict_df.empty:
        # Focus on active uptake reactions that were removed
        active_removed = conflict_df[conflict_df['was_active_uptake']]
        if not active_removed.empty:
            f.write(f"  Active uptake reactions removed: {len(active_removed)}\n\n")
            f.write(f"  {'Sample':<25s} {'MAR_ID':<15s} {'Nutrient':<20s} {'Cause':<18s} {'Scores'}\n")
            f.write("  " + "-" * 110 + "\n")
            for _, row in active_removed.iterrows():
                nutrient = row['nutrient_label'] if row['nutrient_label'] else row['rxn_name']
                scores = row['gene_scores'] if row['gene_scores'] else ''
                f.write(f"  {row['sample']:<25s} {row['rxn_id']:<15s} "
                        f"{str(nutrient)[:19]:<20s} {row['removal_cause']:<18s} {scores}\n")

        # Summary
        cause_counts = conflict_df['removal_cause'].value_counts()
        f.write(f"\n  Removal cause summary:\n")
        for cause, count in cause_counts.items():
            f.write(f"    {cause}: {count}\n")
    else:
        f.write("  No removed exchange reactions to analyze.\n")

    f.write("\n\n" + "=" * 80 + "\n")
    f.write("END OF REPORT\n")
    f.write("=" * 80 + "\n")

print(f"  Report saved: {report_path}")
print("\nDone.")
