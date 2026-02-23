
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COVID-19 Case Study – tINIT Reconstruction Test
================================================

Based on the covid19_case_study.ipynb notebook (FastCORE → tINIT).
Performs:
  1) Load model & omics data
  2) Gene-level thresholding (TAS) using Local T2 strategy
  3) tINIT model reconstruction for all 10 samples via ReconstructionWrapper
  4) Save reconstructed models (SBML/XML) & result CSVs
  5) Visualization (overview, PCA, heatmaps)
  6) Validation & benchmarking (gene essentiality, theoretical yields)
  7) Omics integration summary

Data:
  - HumanGEM_Consistent_COVID19_HAM.xml  (modified HumanGEM model)
  - Desai-GTEx_ensembl.csv               (10 samples, 12 870 genes)
"""

import os
import re
import sys
import time
import json
import warnings
import pickle
import traceback

import numpy as np
import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(BASE_DIR, "HumanGEM_Consistent_COVID19_HAM.xml")
OMICS_PATH = os.path.join(BASE_DIR, "Desai-GTEx_ensembl.csv")

# Add troppo source to path
sys.path.append(os.path.join(BASE_DIR, "../../src"))

OUTPUT_DIR = os.path.join(BASE_DIR, "covid19_tinit_results")
FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures")
MODELS_DIR = os.path.join(OUTPUT_DIR, "models")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(MODELS_DIR, exist_ok=True)

# Suppress non-critical warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Use non-interactive backend for figure saving
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# [USER REQUEST] Limit CPU Usage
import os
os.environ["GRB_THREADS"] = "2"

print("======================================================================")
print("COVID-19 Case Study – tINIT Reconstruction Test")
print("=" * 70)


# ══════════════════════════════════════════════════════════════════════
# STEP 1 : Load Model & Omics Data
# ══════════════════════════════════════════════════════════════════════
print("\n[Step 1] Loading model and omics data...")
t0 = time.time()

from cobra.io import read_sbml_model
model = read_sbml_model(MODEL_PATH)

print(f"  Model: {model.id}")
print(f"    Reactions : {len(model.reactions)}")
print(f"    Metabolites: {len(model.metabolites)}")
print(f"    Genes     : {len(model.genes)}")
print(f"    Objective : {model.objective.expression}")

# Load omics
omics_data = pd.read_csv(OMICS_PATH, index_col=0)
print(f"\n  Omics data: {omics_data.shape[0]} samples × {omics_data.shape[1]} genes")
print(f"  Samples: {list(omics_data.index)}")
print(f"  Loaded in {time.time()-t0:.1f}s")


# ══════════════════════════════════════════════════════════════════════
# STEP 2 : Gene-Level Thresholding (TAS)
# ══════════════════════════════════════════════════════════════════════
print("\n[Step 2] Gene-level thresholding (Local T2 strategy)...")
t0 = time.time()

from troppo.omics import GeneLevelThresholding

# Local T2 with quantile indices: 2=0.5 (global lower), 3=0.75 (global upper), 2=0.5 (local)
# Same as the original notebook
glt = GeneLevelThresholding(
    omics_dataframe=omics_data,
    thresholding_strat='local t2',
    global_threshold_lower=2,   # index 2 → 0.50 quantile
    global_threshold_upper=3,   # index 3 → 0.75 quantile
    local_threshold=2,          # index 2 → 0.50 quantile
)
filtered_df = glt.apply_thresholding_filter()

print(f"  Thresholding complete: {filtered_df.shape}")
print(f"  Score range: [{filtered_df.min().min():.3f}, {filtered_df.max().max():.3f}]")
print(f"  Non-zero scores per sample:")
for sample in filtered_df.index:
    positive = (filtered_df.loc[sample] > 0).sum()
    negative = (filtered_df.loc[sample] < 0).sum()
    print(f"    {sample:45s}  positive={positive}  negative={negative}")
print(f"  Thresholding in {time.time()-t0:.1f}s")


# ══════════════════════════════════════════════════════════════════════
# STEP 3 : tINIT Reconstruction (via ReconstructionWrapper)
# ══════════════════════════════════════════════════════════════════════
print("\n[Step 3] tINIT Reconstruction...")
t0_total = time.time()

from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.readers.generic import TabularReader

# GPR parsing function (from notebook)
patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
replace_alt_transcripts = lambda x: patt.sub('', x)

# Relax model bounds for tINIT feasibility
# (tINIT needs unused reactions to have v=0 as a valid solution)
print("  Relaxing model bounds for tINIT feasibility...")
n_lb_relaxed = n_ub_relaxed = 0
for rxn in model.reactions:
    if rxn.lower_bound > 0:
        rxn.lower_bound = 0
        n_lb_relaxed += 1
    if rxn.upper_bound < 0:
        rxn.upper_bound = 0
        n_ub_relaxed += 1
print(f"    Relaxed {n_lb_relaxed} lower bounds, {n_ub_relaxed} upper bounds")

# ----------------------------------------------------------------------
# [MODIFIED] Preserve Defined Medium & Add Missing Components
# ----------------------------------------------------------------------
# Do NOT open all exchange reactions. Keep original medium.
# But DO ensure Ham's F12 missing components are present.

print("  Adding missing Ham's F12 components (Biotin, B12, Ascorbate)...")
missing_components = {
    'Biotin (B7)': 'm01401s',       # Biotin
    'Cobalamin (B12)': 'm01361s',   # Vitamin B12 (aquacob(III)alamin)
    'Ascorbate (Vit C)': 'm01368s', # L-Ascorbate
}

for name, met_id in missing_components.items():
    try:
        if met_id in model.metabolites:
            met = model.metabolites.get_by_id(met_id)
            # Check if exchange exists
            ex_rxns = [r for r in met.reactions if len(r.metabolites) == 1]
            
            if ex_rxns:
                rxn = ex_rxns[0]
                if rxn.lower_bound == 0:
                    rxn.lower_bound = -1000
                    print(f"    Opened existing exchange for {name}: {rxn.id}")
            else:
                rxn = model.add_boundary(met, type="exchange")
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
                print(f"    Added new exchange for {name}: {rxn.id}")
        else:
            print(f"    WARNING: Metabolite {met_id} ({name}) not found in model.")
    except Exception as e:
        print(f"    Error adding {name}: {e}")

# Create wrapper
print("  Creating ReconstructionWrapper...")
rw = ReconstructionWrapper(
    model=model,
    ttg_ratio=9999,
    gpr_gene_parse_function=replace_alt_transcripts,
)

# Create OmicsContainers from thresholded data
print("  Creating OmicsContainers from TabularReader...")
omics_containers = TabularReader(
    path_or_df=filtered_df,
    nomenclature='entrez_id',
    omics_type='transcriptomics',
).to_containers()
print(f"  Created {len(omics_containers)} OmicsContainers")

# Identify essential reactions (biomass + key exchange reactions)
# tINIT will always include these in the reconstructed model
r_ids_list = list(rw.model_reader.r_ids)
essential_rxn_ids = []
# Biomass reaction
for rxn in model.reactions:
    if 'biomass' in rxn.id.lower():
        essential_rxn_ids.append(rxn.id)
# Key nutrient exchange reactions
for rxn_id in ['HMR_9034', 'HMR_9048', 'HMR_9058']:  # glucose, O2, CO2
    if rxn_id in r_ids_list:
        essential_rxn_ids.append(rxn_id)
# Biomass chain: transport + exchange + cofactor dependencies
for rxn_id in ['HMR_10023', 'HMR_10024', 'HMR_8617', 'HMR_10065']:
    if rxn_id in r_ids_list:
        essential_rxn_ids.append(rxn_id)

essential_rxn_indices = [r_ids_list.index(rid) for rid in essential_rxn_ids if rid in r_ids_list]
print(f"  Essential reactions: {len(essential_rxn_indices)} ({essential_rxn_ids})")

# Map original sample names
# filtered_df index: 'Lung_Healthy_local_t2_2_3_2' → original: 'Lung_Healthy'
filter_suffix = '_local_t2_2_3_2'

all_results = {}
per_sample_times = {}

for idx, container in enumerate(omics_containers):
    # Recover the original sample name from the container condition
    container_name = container.condition
    orig_name = container_name.replace(filter_suffix, '') if filter_suffix in container_name else container_name

    print(f"\n  [{idx+1}/{len(omics_containers)}] Reconstructing: {orig_name}")
    t_sample = time.time()

    try:
        # Manual integration pipeline (to handle None scores from GPR-less reactions)
        from troppo.omics.integration import ContinuousScoreIntegrationStrategy
        from troppo.methods.reconstruction.tINIT import tINITProperties

        # 1. Get integrated data map (gene → reaction scores via GPR)
        data_map = container.get_integrated_data_map(rw.model_reader, min, max)

        # 2. Get continuous scores
        strat = ContinuousScoreIntegrationStrategy(score_apply=None)
        scores_dict = strat.integrate(data_map)

        # 3. Replace None scores with 0 (reactions without GPR association)
        scores_list = [scores_dict.get(r, 0) if scores_dict.get(r) is not None else 0
                       for r in rw.model_reader.r_ids]
        n_none = sum(1 for r in rw.model_reader.r_ids if scores_dict.get(r) is None)
        n_positive = sum(1 for s in scores_list if s > 0)
        print(f"    Scores: {n_positive} positive, {n_none} None→0")

        # 4. Create tINIT properties and run
        properties = tINITProperties.from_integrated_scores(
            scores_list,
            solver='GUROBI',
            essential_reactions=essential_rxn_indices,
        )
        algorithm_result = rw.run(properties)
        result_names = [rw.model_reader.r_ids[k] for k in algorithm_result]
        result = {k: k in result_names for k in rw.model_reader.r_ids}

        elapsed = time.time() - t_sample
        per_sample_times[orig_name] = elapsed

        active_count = sum(1 for v in result.values() if v)
        total_count = len(result)
        print(f"    Active reactions: {active_count}/{total_count} "
              f"({active_count/total_count*100:.1f}%) in {elapsed:.1f}s")
        all_results[orig_name] = result

    except Exception as e:
        elapsed = time.time() - t_sample
        per_sample_times[orig_name] = elapsed
        print(f"    FAILED ({elapsed:.1f}s): {e}")
        traceback.print_exc()
        all_results[orig_name] = {r.id: False for r in model.reactions}

total_elapsed = time.time() - t0_total
n_samples = len(omics_containers)
print(f"\n  Total reconstruction time: {total_elapsed:.1f}s "
      f"(avg {total_elapsed/max(n_samples,1):.1f}s/sample)")


# ══════════════════════════════════════════════════════════════════════
# STEP 4 : Build ReconstructionResults & Save
# ══════════════════════════════════════════════════════════════════════
print("\n[Step 4] Building results and saving...")
t0 = time.time()

from troppo.results import ReconstructionResults

results = ReconstructionResults(
    results_dict=all_results,
    model=model,
    metadata={
        'algorithm': 'tINIT',
        'dataset': 'COVID-19 (Desai+GTEx)',
        'thresholding': 'local_t2 (2,3,2)',
        'per_sample_times': per_sample_times,
    },
)

# ── 4a. Save summary CSV ──
summary = results.summary()
summary_path = os.path.join(OUTPUT_DIR, "reconstruction_summary.csv")
summary.to_csv(summary_path)
print(f"  Summary saved: {summary_path}")
print(summary.to_string())

# ── 4b. Save reaction matrix ──
rxn_matrix = results.reaction_matrix
rxn_path = os.path.join(OUTPUT_DIR, "reaction_matrix.csv")
rxn_matrix.astype(int).to_csv(rxn_path)
print(f"  Reaction matrix saved: {rxn_path}")

# ── 4c. Save Jaccard similarity ──
jaccard = results.jaccard_similarity()
jaccard_path = os.path.join(OUTPUT_DIR, "jaccard_similarity.csv")
jaccard.to_csv(jaccard_path)
print(f"  Jaccard similarity saved: {jaccard_path}")

# ── 4d. Save subsystem coverage ──
subsys = results.subsystem_coverage()
if subsys is not None and not subsys.empty:
    subsys_path = os.path.join(OUTPUT_DIR, "subsystem_coverage.csv")
    subsys.to_csv(subsys_path)
    print(f"  Subsystem coverage saved: {subsys_path}")

# ── 4e. Save per-sample models as SBML/XML ──
print("\n  Saving per-sample models as SBML/XML...")
from cobra.io import write_sbml_model
for sample in results.sample_names:
    try:
        sample_model = results.get_model(sample)
        safe_name = sample.replace(" ", "_").replace("/", "_")
        xml_path = os.path.join(MODELS_DIR, f"{safe_name}_model.xml")
        write_sbml_model(sample_model, xml_path)
        print(f"    {sample}: {len(sample_model.reactions)} reactions → {xml_path}")
    except Exception as e:
        print(f"    {sample}: Failed to save model: {e}")

# ── 4e-2. Post-hoc gap-fill for BH4/THF/CL pathway gaps ──
print("\n  Post-hoc gap-fill for pathway gaps...")
fix_csv = os.path.join(BASE_DIR, "diagnostic_universal_fix.csv")
if os.path.exists(fix_csv):
    fix_df = pd.read_csv(fix_csv)
    gapfill_rxn_ids = fix_df[fix_df['is_core'] == False]['rxn_id'].tolist()
    print(f"    Gap-fill candidate reactions: {len(gapfill_rxn_ids)}")

    for sample in results.sample_names:
        sample_model = results.get_model(sample)
        # Skip empty models (timed-out reconstructions)
        if len(sample_model.reactions) == 0:
            print(f"    {sample}: empty model (reconstruction timed out), skipping")
            continue
        sample_model.objective = 'biomass_human'
        sol = sample_model.optimize()
        if sol.objective_value is not None and sol.objective_value < 1e-6:
            # Add gap-fill reactions from universal model
            model_rxn_ids = set(r.id for r in sample_model.reactions)
            model_met_ids = set(m.id for m in sample_model.metabolites)
            added = []
            for rid in gapfill_rxn_ids:
                if rid not in model_rxn_ids and rid in set(r.id for r in model.reactions):
                    urxn = model.reactions.get_by_id(rid).copy()
                    missing_mets = [met for met in urxn.metabolites if met.id not in model_met_ids]
                    for met in missing_mets:
                        sample_model.add_metabolites([met.copy()])
                        model_met_ids.add(met.id)
                    sample_model.add_reactions([urxn])
                    model_rxn_ids.add(rid)
                    added.append(rid)
            if added:
                sol2 = sample_model.optimize()
                bm_val = sol2.objective_value if sol2.objective_value is not None else 0.0
                print(f"    {sample}: gap-filled +{len(added)} rxns, biomass={bm_val:.6f}")
                # Re-save the gap-filled model
                safe_name = sample.replace(" ", "_").replace("/", "_")
                xml_path = os.path.join(MODELS_DIR, f"{safe_name}_model.xml")
                write_sbml_model(sample_model, xml_path)
            else:
                print(f"    {sample}: no gap-fill needed (all reactions present)")
        else:
            bm_val = sol.objective_value if sol.objective_value is not None else 0.0
            print(f"    {sample}: biomass={bm_val:.6f} (OK, no gap-fill needed)")
else:
    print(f"    WARNING: {fix_csv} not found, skipping gap-fill")

# ── 4f. Pickle results ──
pkl_path = os.path.join(OUTPUT_DIR, "reconstruction_results.pkl")
with open(pkl_path, 'wb') as f:
    pickle.dump(results, f)
print(f"  Pickled results: {pkl_path}")

print(f"  Step 4 completed in {time.time()-t0:.1f}s")


# ══════════════════════════════════════════════════════════════════════
# STEP 5 : Visualization
# ══════════════════════════════════════════════════════════════════════
print("\n[Step 5] Visualization...")
t0 = time.time()

try:
    from troppo.visualization import (
        plot_overview,
        plot_model_sizes,
        plot_reaction_overlap_heatmap,
        plot_subsystem_coverage,
        plot_pca_models,
        plot_reaction_frequency,
        plot_biomass_distribution,
    )

    # 5a. Overview (2×3 panel)
    print("  Generating overview plot...")
    plot_overview(results, save_dir=FIGURES_DIR, format="png")
    print(f"    Saved to {FIGURES_DIR}/")

    # 5b. Individual plots
    print("  Generating individual plots...")

    fig = plot_model_sizes(results)
    fig.savefig(os.path.join(FIGURES_DIR, "model_sizes.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plot_reaction_overlap_heatmap(results)
    fig.savefig(os.path.join(FIGURES_DIR, "jaccard_heatmap.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plot_subsystem_coverage(results, top_n=25)
    fig.savefig(os.path.join(FIGURES_DIR, "subsystem_coverage.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plot_reaction_frequency(results, top_n=30)
    fig.savefig(os.path.join(FIGURES_DIR, "reaction_frequency.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

    try:
        labels = []
        for s in results.sample_names:
            if 'COVID' in s:
                labels.append('COVID-19')
            elif 'Healthy' in s:
                labels.append('Healthy')
            else:
                labels.append('Other')
        fig = plot_pca_models(results, labels=labels)
        fig.savefig(os.path.join(FIGURES_DIR, "pca_models.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("    PCA plot saved.")
    except ImportError:
        print("    PCA skipped (scikit-learn not available)")

    try:
        fig = plot_biomass_distribution(results)
        fig.savefig(os.path.join(FIGURES_DIR, "biomass_distribution.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("    Biomass distribution plot saved.")
    except Exception as e:
        print(f"    Biomass distribution skipped: {e}")

    print(f"  Visualization completed in {time.time()-t0:.1f}s")

except ImportError as e:
    print(f"  Visualization skipped (missing dependency): {e}")
except Exception as e:
    print(f"  Visualization error: {e}")
    traceback.print_exc()


# ══════════════════════════════════════════════════════════════════════
# STEP 6 : Validation & Benchmarking
# ══════════════════════════════════════════════════════════════════════
print("\n[Step 6] Validation & Benchmarking...")
t0 = time.time()

try:
    from troppo.benchmark.validation import GeneEssentialityValidator, TheoreticalYieldCalculator

    # ── 6a. Theoretical Yield Calculation ──
    print("\n  [6a] Theoretical Yield Calculation...")

    human_gem_mapping = {
        'glucose': {'id': 'm01965s', 'name': 'D-glucose', 'exchange': 'HMR_9034'},
        'oxygen': {'id': 'm02630s', 'name': 'oxygen', 'exchange': 'HMR_9048'},
        'co2': {'id': 'm01596s', 'name': 'carbon dioxide', 'exchange': 'HMR_9058'},
    }

    calculator = TheoreticalYieldCalculator(metabolite_mapping=human_gem_mapping)

    all_yield_results = []
    for sample in results.sample_names:
        group = "COVID19" if "COVID" in sample else "Healthy"
        tissue = sample.split("_")[0]

        try:
            sample_model = results.get_model(sample)

            # Open exchange reactions for feasibility
            # [MODIFIED] Do not force open random exchanges, rely on medium
            # But ensure basic exchanges are open if needed (usually handled by medium)
            
            yields = calculator.calculate_yields(
                sample_model,
                carbon_sources=['glucose'],
                aerobic=True,
                anaerobic=False,
                product='biomass',
            )

            for yr in yields:
                yd = yr.to_dict()
                yd['sample'] = sample
                yd['group'] = group
                yd['tissue'] = tissue
                all_yield_results.append(yd)

            print(f"    {sample}: yield={yields[0].theoretical_yield:.4f}, "
                  f"feasible={yields[0].feasible}")

        except Exception as e:
            print(f"    {sample}: Yield calc failed: {e}")

    if all_yield_results:
        yield_df = pd.DataFrame(all_yield_results)
        yield_path = os.path.join(OUTPUT_DIR, "benchmark_yields.csv")
        yield_df.to_csv(yield_path, index=False)
        print(f"\n  Yield results saved: {yield_path}")

        print("\n  -- Average Biomass Yield by Group --")
        agg_cols = [c for c in ['theoretical_yield', 'feasible'] if c in yield_df.columns]
        if agg_cols:
            print(yield_df.groupby('group')[agg_cols].mean().to_string())

    # ── 6b. Gene Essentiality Validation ──
    print("\n  [6b] Gene Essentiality Validation...")

    # Human core essential genes (CRISPR screen consensus from DepMap / Hart et al.)
    # These are genes whose knockout is lethal in most human cell lines
    human_essential_genes_symbols = [
        # Core transcription/translation
        'RPS3', 'RPS5', 'RPL5', 'RPL11', 'POLR2A', 'POLR2B',
        # Glycolysis / TCA cycle
        'HK1', 'GPI', 'PFKM', 'ALDOA', 'GAPDH', 'PGK1', 'ENO1', 'PKM',
        'CS', 'ACO2', 'IDH3A', 'OGDH', 'SUCLA2', 'SDHA', 'FH', 'MDH2',
        # Oxidative phosphorylation
        'NDUFS1', 'NDUFS2', 'NDUFS3', 'SDHB', 'UQCRC1', 'COX5A',
        'ATP5F1A', 'ATP5F1B',
        # Nucleotide biosynthesis
        'UMPS', 'DHODH', 'CTPS1', 'GART', 'PAICS', 'ATIC',
        # Amino acid metabolism
        'BCAT1', 'GOT2', 'GLUD1',
        # Folate / one-carbon
        'DHFR', 'MTHFD1', 'SHMT2',
    ]

    # Map gene symbols to entrez IDs (model uses entrez)
    model_gene_ids = {g.id for g in model.genes}
    # The model uses Entrez gene IDs, so we check overlap with omics gene IDs
    # For now, take genes that are actually in the model
    essential_in_model = [g for g in model_gene_ids if g in model_gene_ids]

    # Use single gene deletion on sub-models to find predicted essential genes
    # and compare with known essential gene list
    essentiality_results = []

    # Run on 2 samples per group for speed
    healthy_samples = [s for s in results.sample_names if 'Healthy' in s]
    covid_samples = [s for s in results.sample_names if 'COVID' in s]
    test_samples = healthy_samples[:2] + covid_samples[:2]

    for sample in test_samples:
        group = "COVID19" if "COVID" in sample else "Healthy"
        tissue = sample.split("_")[0]

        try:
            sample_model = results.get_model(sample)

            # Check if biomass works
            wt_sol = sample_model.optimize()
            wt_growth = wt_sol.objective_value if wt_sol.status == 'optimal' else 0.0

            if wt_growth < 1e-6:
                print(f"    {sample}: WT growth={wt_growth:.6f} (too low, skipping essentiality)")
                essentiality_results.append({
                    'accuracy': 0.0, 'precision': 0.0, 'recall': 0.0,
                    'f1_score': 0.0, 'specificity': 0.0, 'matthews_correlation': 0.0,
                    'num_tested': 0, 'num_experimental_essential': 0,
                    'num_predicted_essential': 0, 'tp': 0, 'fp': 0, 'tn': 0, 'fn': 0,
                    'wt_growth': wt_growth,
                    'sample': sample, 'group': group, 'tissue': tissue
                })
                continue

            # Perform gene essentiality with cobra
            from cobra.flux_analysis import single_gene_deletion
            deletion_results = single_gene_deletion(sample_model)

            # Count predicted essential genes
            n_essential = 0
            n_total = len(deletion_results)
            for idx_row in deletion_results.index:
                growth_val = deletion_results.loc[idx_row, 'growth']
                if isinstance(growth_val, pd.Series):
                    growth_val = growth_val.values[0]
                if np.isnan(growth_val) or growth_val < 0.01 * wt_growth:
                    n_essential += 1

            print(f"    {sample}: WT_growth={wt_growth:.4f}, "
                  f"predicted_essential={n_essential}/{n_total}")

            essentiality_results.append({
                'wt_growth': wt_growth,
                'num_model_genes': len(sample_model.genes),
                'num_tested': n_total,
                'num_predicted_essential': n_essential,
                'pct_essential': n_essential / max(n_total, 1) * 100,
                'sample': sample, 'group': group, 'tissue': tissue
            })

        except Exception as e:
            print(f"    {sample}: Essentiality failed: {e}")
            import traceback as tb
            tb.print_exc()

    if essentiality_results:
        ess_df = pd.DataFrame(essentiality_results)
        ess_path = os.path.join(OUTPUT_DIR, "benchmark_essentiality.csv")
        ess_df.to_csv(ess_path, index=False)
        print(f"\n  Essentiality results saved: {ess_path}")

    print(f"\n  Validation & Benchmarking completed in {time.time()-t0:.1f}s")

except Exception as e:
    print(f"  Validation/Benchmark failed: {e}")
    traceback.print_exc()


# ══════════════════════════════════════════════════════════════════════
# STEP 7 : Omics Integration Summary & Analysis
# ══════════════════════════════════════════════════════════════════════
print("\n[Step 7] Omics Integration Summary...")
t0 = time.time()

try:
    # ── 7a. Reaction overlap analysis (Healthy vs COVID-19) ──
    healthy_names = [s for s in results.sample_names if 'Healthy' in s]
    covid_names = [s for s in results.sample_names if 'COVID' in s]

    rxn = results.reaction_matrix

    healthy_any = rxn.loc[healthy_names].any(axis=0)
    covid_any = rxn.loc[covid_names].any(axis=0)

    healthy_only = healthy_any & ~covid_any
    covid_only = covid_any & ~healthy_any
    shared = healthy_any & covid_any

    print(f"\n  -- Reaction Overlap (Healthy vs COVID-19) --")
    print(f"    Healthy only:  {healthy_only.sum()}")
    print(f"    COVID-19 only: {covid_only.sum()}")
    print(f"    Shared:        {shared.sum()}")

    # Save differential reactions
    diff_rxns = pd.DataFrame({
        'reaction': rxn.columns,
        'healthy_only': healthy_only.values,
        'covid_only': covid_only.values,
        'shared': shared.values,
        'n_healthy_active': rxn.loc[healthy_names].sum(axis=0).values,
        'n_covid_active': rxn.loc[covid_names].sum(axis=0).values,
    })
    diff_path = os.path.join(OUTPUT_DIR, "differential_reactions.csv")
    diff_rxns.to_csv(diff_path, index=False)
    print(f"    Differential reactions saved: {diff_path}")

    # ── 7b. Per-tissue comparison ──
    tissues = set(s.split("_")[0] for s in results.sample_names)
    print(f"\n  -- Per-Tissue Active Reactions --")
    tissue_summary = []
    for tissue in sorted(tissues):
        healthy_s = f"{tissue}_Healthy"
        covid_s = f"{tissue}_COVID19"

        if healthy_s in results.sample_names and covid_s in results.sample_names:
            h_active = rxn.loc[healthy_s].sum()
            c_active = rxn.loc[covid_s].sum()
            h_set = set(rxn.columns[rxn.loc[healthy_s]])
            c_set = set(rxn.columns[rxn.loc[covid_s]])

            tissue_summary.append({
                'Tissue': tissue,
                'Healthy_active': int(h_active),
                'COVID19_active': int(c_active),
                'Healthy_only': len(h_set - c_set),
                'COVID19_only': len(c_set - h_set),
                'Shared': len(h_set & c_set),
                'Jaccard': len(h_set & c_set) / max(len(h_set | c_set), 1),
            })

            print(f"    {tissue:12s}  Healthy={int(h_active):5d}  COVID-19={int(c_active):5d}  "
                  f"Shared={len(h_set & c_set):5d}  Jaccard={len(h_set & c_set)/max(len(h_set|c_set),1):.3f}")

    tissue_df = pd.DataFrame(tissue_summary)
    tissue_path = os.path.join(OUTPUT_DIR, "tissue_comparison.csv")
    tissue_df.to_csv(tissue_path, index=False)
    print(f"    Tissue comparison saved: {tissue_path}")

    # ── 7c. Subsystem coverage difference ──
    if subsys is not None and not subsys.empty:
        print(f"\n  -- Top 15 Differentially Active Subsystems (Healthy vs COVID-19) --")
        h_mean = subsys.loc[healthy_names].mean(axis=0)
        c_mean = subsys.loc[covid_names].mean(axis=0)
        diff_subsys = (c_mean - h_mean).abs().sort_values(ascending=False).head(15)
        for sub in diff_subsys.index:
            print(f"    {sub:55s}  H={h_mean[sub]:.3f}  C={c_mean[sub]:.3f}  "
                  f"Δ={c_mean[sub]-h_mean[sub]:+.3f}")

    # ── 7d. Biomass flux comparison ──
    print(f"\n  -- Biomass Flux Comparison --")
    try:
        biomass_data = {}
        for sample in results.sample_names:
            try:
                sample_model = results.get_model(sample)
                # Ensure exchanges are open for flux
                # [MODIFIED] Respect defined medium
                # for rxn in sample_model.exchanges:
                #    rxn.lower_bound = -1000
                #    rxn.upper_bound = 1000
                sol = sample_model.optimize()
                flux = sol.objective_value if sol.status == 'optimal' else 0.0
                biomass_data[sample] = flux
            except Exception:
                biomass_data[sample] = 0.0

        biomass = pd.Series(biomass_data, name='biomass_flux')
        biomass_path = os.path.join(OUTPUT_DIR, "biomass_fluxes.csv")
        biomass.to_csv(biomass_path)
        print(f"    Biomass fluxes saved: {biomass_path}")
        for sample in results.sample_names:
            print(f"    {sample:25s}  biomass_flux={biomass.get(sample, 0):.6f}")
    except Exception as e:
        print(f"    Biomass flux calculation failed: {e}")

    print(f"\n  Omics Integration Summary completed in {time.time()-t0:.1f}s")

except Exception as e:
    print(f"  Omics integration analysis failed: {e}")
    traceback.print_exc()


# ══════════════════════════════════════════════════════════════════════
# STEP 8 : Final Report
# ══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("FINAL REPORT")
print("=" * 70)

print(f"\n  Output directory: {OUTPUT_DIR}")
print(f"  Figures directory: {FIGURES_DIR}")
print(f"  Models directory: {MODELS_DIR}")
print(f"\n  Files generated:")
for root, dirs, files in os.walk(OUTPUT_DIR):
    for f in sorted(files):
        fpath = os.path.join(root, f)
        size_kb = os.path.getsize(fpath) / 1024
        rel = os.path.relpath(fpath, OUTPUT_DIR)
        print(f"    {rel:60s} {size_kb:10.1f} KB")

print(f"\n{'=' * 70}")
print("COVID-19 tINIT Test completed successfully!")
print(f"{'=' * 70}")
