"""
MCF7 Context-Specific GEM Reconstruction via tINIT
===================================================

Workflow
--------
  Stage 1 — Reference model essential task evaluation
  Stage 2 — FVA-based biomass-essential reaction detection
  Stage 3 — tINIT reconstruction with task validation + essential reactions
  Stage 4 — Model visualization (size, subsystem coverage, reaction frequency)
  Stage 5 — Biomass flux via FBA on the context-specific model
  Stage 6 — Theoretical yield calculation (glucose & glutamine, aerobic/anaerobic)

Scientific background (Robinson et al. 2020)
-------------------------------------------
  57 essential metabolic tasks span five categories:
    ER  Energy & Redox        — ATP/GTP, TCA, OXPHOS
    BS  Biosynthesis          — nucleotides, cholesterol, phospholipids, vitamins
    SU  Substrate Utilization — amino acids, fatty acids, choline, inositol
    IC  Internal Conversions  — glycolysis, acetyl-CoA, TCA intermediates
    GR  Growth                — biomass production on Ham's F12 medium

  The GR1 task trivially passes via reverse flux of HMR_10024 (bidirectional
  biomass exchange).  FVA-based auto-detection (derive_biomass_essential_reactions)
  identifies the minimum set of reactions required for genuine biomass production
  and forces them into the tINIT reconstruction via essential_reactions.

Usage
-----
  python test_tinit_task_validation.py [--sample SAMPLE] [--no-plots]

References
----------
  Robinson JL et al. (2020) An atlas of human metabolism.
    Science Signaling 13(624):eaaz1482.
    doi: 10.1126/scisignal.aaz1482
"""

import os
import sys
import re
import warnings
import argparse
from collections import defaultdict

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

import cobra
import numpy as np
import pandas as pd

from troppo import TINITReconstructor
from troppo.tasks.core import TaskEvaluator
from troppo.benchmark.validation import TheoreticalYieldCalculator

# ── Path configuration ─────────────────────────────────────────────────────────
DATA_DIR    = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
MODEL_PATH  = os.path.join(DATA_DIR, "Human-GEM_consistent.xml")
EXPR_PATH   = os.path.join(DATA_DIR, "CCLE_expression_ACH-000019_scores.csv")
TASK_PATH   = os.path.join(DATA_DIR, "metabolicTasks_Essential.xlsx")
BIOMASS_ID  = "biomass_human"

CATEGORY_DESCRIPTIONS = {
    "ER": "Energy & Redox (ATP rephosphorylation, TCA, OXPHOS)",
    "BS": "Biosynthesis (nucleotides, cholesterol, phospholipids, vitamins)",
    "SU": "Substrate Utilization (amino acids, fatty acids, choline, inositol)",
    "IC": "Internal Conversions (glycolysis, acetyl-CoA, TCA intermediates)",
    "GR": "Growth (biomass production on Ham's F12 medium)",
}

SEP  = "=" * 64
DASH = "─" * 60


# ── CCLE expression data loader ────────────────────────────────────────────────
def load_ccle(path: str) -> pd.DataFrame:
    """Load CCLE score CSV as a genes × samples DataFrame (global configs only)."""
    df = pd.read_csv(path)
    gene_cols = [c for c in df.columns if "ENSG" in str(c)]
    meta_cols = [c for c in df.columns if c.startswith("Unnamed")]
    df_expr   = df[gene_cols].copy()
    df_expr.index = df[meta_cols].apply(
        lambda r: "_".join(str(v) for v in r), axis=1
    )
    df_T = df_expr.T
    df_T.index = [
        re.search(r"(ENSG\d+)", str(c)).group(1)
        if re.search(r"(ENSG\d+)", str(c)) else c
        for c in df_T.index
    ]
    return df_T[[c for c in df_T.columns if c.startswith("global_")]]


def get_category(name: str) -> str:
    """Map task name to category code (e.g. 'ER1' -> 'ER')."""
    for cat in CATEGORY_DESCRIPTIONS:
        if name.startswith(cat):
            return cat
    return "OTHER"


# ── Stage 1: Reference model task evaluation ──────────────────────────────────
def stage1_reference(recon, tasks: list) -> dict:
    """Evaluate essential tasks on the fully active reference model."""
    print(f"\n{SEP}")
    print("  STAGE 1: Reference Model Task Evaluation")
    print(f"  Human-GEM_consistent — {len(tasks)} compatible essential tasks")
    print(SEP)

    tev = TaskEvaluator(model=recon._model, tasks=tasks, solver="GUROBI")
    raw = tev.batch_evaluate(bound_changes=[{}], threads=1)

    results   = {}
    cat_stats = defaultdict(lambda: [0, 0])  # [passed, total]

    for task in tasks:
        truth = raw[0].get(task.name, (False, {}, None))[0]
        results[task.name] = truth
        cat = get_category(task.name)
        cat_stats[cat][1] += 1
        if truth:
            cat_stats[cat][0] += 1

    for cat in ["ER", "BS", "SU", "IC", "GR"]:
        p, t = cat_stats[cat]
        if t == 0:
            continue
        bar = "#" * p + "-" * (t - p)
        print(f"  {cat}  [{bar}]  {p}/{t}  {CATEGORY_DESCRIPTIONS[cat]}")

    passed = sum(results.values())
    total  = len(results)
    print(f"  {DASH}")
    print(f"  TOTAL: {passed}/{total} passed")
    return results


# ── Stage 2: FVA-based biomass essential reaction detection ───────────────────
def stage2_fva_essential(recon) -> tuple:
    """Identify biomass-essential reactions via FVA on the reference model.

    Uses recon.derive_biomass_essential_reactions() (see README.rst for
    the three-step pattern: detect -> configure -> run).
    """
    print(f"\n{SEP}")
    print("  STAGE 2: FVA-based Biomass Essential Reaction Detection")
    print(f"  Biomass reaction: {BIOMASS_ID}  |  min_fraction = 1%")
    print(SEP)

    essential, fva_df = recon.derive_biomass_essential_reactions(
        biomass_rxn_id=BIOMASS_ID,
        min_fraction=0.01,
    )
    print(f"  FVA essential reactions identified: {len(essential)}")
    return essential, fva_df


# ── Stage 3: tINIT reconstruction + task validation ──────────────────────────
def stage3_reconstruction(recon, sample: str, essential: list) -> object:
    """Run tINIT with FVA-derived essential reactions and validate essential tasks.

    Three-step pattern (from README.rst Biomass Flux section):
      1. derive_biomass_essential_reactions() — already done in stage2
      2. configure(essential_reactions=essential)
      3. run(validate_tasks=True, gap_fill=True)
    """
    print(f"\n{SEP}")
    print(f"  STAGE 3: tINIT Reconstruction  (sample: {sample})")
    print(f"  essential_reactions: {len(essential)}  |  gap_fill: True  |  validate_tasks: True")
    print(SEP)

    recon.configure(essential_reactions=essential)
    results_obj = recon.run(
        samples=[sample],
        gap_fill=True,
        validate_tasks=True,
        task_file=TASK_PATH,
        stream_tasks=False,
    )

    rm       = results_obj.reaction_matrix
    n_active = int(rm.loc[sample].sum())
    n_total  = rm.shape[1]
    print(f"  Active reactions : {n_active}/{n_total} ({n_active / n_total * 100:.1f}%)")

    # Confirm essential reactions are present in the reconstruction
    rd      = rm.loc[sample].to_dict()
    missing = [r for r in essential if not rd.get(r, False)]
    if missing:
        print(f"  [!] Essential reactions not active: {len(missing)}")
    else:
        print(f"  [OK] All {len(essential)} essential reactions are active")

    # Task validation results from run()
    task_res = results_obj._task_results.get(sample)
    if task_res:
        t_pass  = sum(1 for v in task_res.results.values() if v)
        t_total = len(task_res.results)
        failed  = [n for n, v in task_res.results.items() if not v]
        print(f"\n  Essential tasks  : {t_pass}/{t_total} passed")
        if failed:
            print(f"  Failed tasks     : {failed}")

        cat_stats = defaultdict(lambda: [0, 0])
        for name, truth in task_res.results.items():
            cat = get_category(name)
            cat_stats[cat][1] += 1
            if truth:
                cat_stats[cat][0] += 1

        print(f"\n  {'Cat':<5} {'Passed':>10}  Description")
        print(f"  {DASH}")
        for cat in ["ER", "BS", "SU", "IC", "GR"]:
            p, t = cat_stats[cat]
            if t == 0:
                continue
            desc = CATEGORY_DESCRIPTIONS.get(cat, "")[:40]
            print(f"  {cat:<5} {p:>3}/{t:<6}  {desc}")

    return results_obj


# ── Stage 4: Model visualization ─────────────────────────────────────────────
def stage4_visualize(results_obj, save_dir: str):
    """Generate and save model visualization plots using troppo.visualization."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from troppo.visualization import (
            plot_model_sizes,
            plot_subsystem_coverage,
            plot_reaction_frequency,
        )
    except ImportError as e:
        print(f"  [!] Visualization skipped (import error): {e}")
        return

    print(f"\n{SEP}")
    print("  STAGE 4: Model Visualization")
    print(SEP)

    os.makedirs(save_dir, exist_ok=True)

    plot_specs = [
        ("model_sizes",        plot_model_sizes,        "model_sizes.png"),
        ("subsystem_coverage", plot_subsystem_coverage,  "subsystem_coverage.png"),
        ("reaction_frequency", plot_reaction_frequency,  "reaction_frequency.png"),
    ]

    for label, fn, fname in plot_specs:
        try:
            out = os.path.join(save_dir, fname)
            fig = fn(results_obj, save_path=out)
            plt.close(fig)
            print(f"  [OK] {label:<22} -> {fname}")
        except Exception as e:
            print(f"  [!!] {label:<22} -> {e}")


# ── Stage 5: Biomass flux ─────────────────────────────────────────────────────
def stage5_biomass(recon, results_obj, sample: str) -> tuple:
    """Compute FBA biomass flux for both the reference and context-specific models."""
    print(f"\n{SEP}")
    print("  STAGE 5: Biomass Flux (FBA)")
    print(SEP)

    # Reference model (full, all reactions active)
    with recon._model:
        recon._model.objective = BIOMASS_ID
        sol_ref = recon._model.optimize()
    v_ref = sol_ref.objective_value if sol_ref.status == "optimal" else 0.0

    # Context-specific model (inactive reactions removed via get_model)
    v_tinit = 0.0
    try:
        ctx = results_obj.get_model(sample)
        ctx.objective = BIOMASS_ID
        sol2  = ctx.optimize()
        v_tinit = sol2.objective_value if sol2.status == "optimal" else 0.0
    except Exception as e:
        print(f"  [!] get_model failed ({e}), using bound-change fallback")
        rd  = results_obj.reaction_matrix.loc[sample].to_dict()
        bc  = {r.id: (0.0, 0.0) for r in recon._model.reactions if not rd.get(r.id, False)}
        with recon._model:
            for rid, bnd in bc.items():
                recon._model.reactions.get_by_id(rid).bounds = bnd
            recon._model.objective = BIOMASS_ID
            sol2 = recon._model.optimize()
        v_tinit = sol2.objective_value if sol2.status == "optimal" else 0.0

    ratio = v_tinit / v_ref * 100 if v_ref > 1e-9 else 0.0
    print(f"  {'Model':<36} {'Biomass (mmol/gDW/h)':>20}  {'vs Ref':>7}")
    print(f"  {DASH}")
    print(f"  {'Reference (Human-GEM_consistent)':<36} {v_ref:>20.4f}   {'100.0%':>7}")
    print(f"  {'tINIT (' + sample + ')':<36} {v_tinit:>20.4f}   {ratio:>6.1f}%")

    return v_ref, v_tinit


# ── Human-GEM exchange reaction discovery ─────────────────────────────────────
def find_exchange_by_name(model: cobra.Model, keywords: list):
    """Return the first exchange reaction whose metabolite name matches any keyword."""
    for rxn in model.exchanges:
        for met in rxn.metabolites:
            for kw in keywords:
                if kw.lower() in met.name.lower():
                    return rxn
    return None


def build_humangem_metabolite_mapping(model: cobra.Model) -> dict:
    """Build a TheoreticalYieldCalculator metabolite mapping for Human-GEM.

    Human-GEM uses HMDB IDs rather than BiGG notation, so the default
    metabolite mapping in TheoreticalYieldCalculator does not apply.
    This function searches exchange reactions by metabolite name keywords
    and returns a mapping dict compatible with TheoreticalYieldCalculator.
    """
    search_targets = {
        "glucose":   ["d-glucose", "glucose"],
        "glutamine": ["l-glutamine", "glutamine"],
        "pyruvate":  ["pyruvate"],
        "lactate":   ["l-lactate", "lactate"],
        "oxygen":    ["oxygen", "o2"],
        "co2":       ["carbon dioxide", "co2"],
    }

    mapping = {}
    for key, keywords in search_targets.items():
        rxn = find_exchange_by_name(model, keywords)
        if rxn is not None:
            met = next(iter(rxn.metabolites))
            mapping[key] = {
                "id":       met.id,
                "name":     met.name,
                "exchange": rxn.id,
            }
    return mapping


# ── Stage 6: Theoretical yield calculation ────────────────────────────────────
def stage6_yields(recon, results_obj, sample: str, save_dir: str) -> pd.DataFrame:
    """Compute maximum theoretical biomass yields for representative carbon sources.

    Follows the benchmark example pattern (benchmark_with_validation_example.py):
      TheoreticalYieldCalculator.calculate_yields() with aerobic=True, anaerobic=True.

    For Human-GEM, a custom metabolite mapping is built from exchange reaction
    names because the model uses HMDB IDs rather than BiGG notation.
    """
    print(f"\n{SEP}")
    print("  STAGE 6: Theoretical Yield Calculation")
    print(SEP)

    # Discover Human-GEM exchange reactions by metabolite name
    mapping = build_humangem_metabolite_mapping(recon._model)
    carbon_sources = [c for c in ["glucose", "glutamine", "pyruvate", "lactate"]
                      if c in mapping]

    if not carbon_sources:
        print("  [!] No carbon source exchange reactions found — skipping.")
        return pd.DataFrame()

    print("  Carbon sources identified in model:")
    for cs in carbon_sources:
        print(f"    {cs:<12}  exchange: {mapping[cs]['exchange']}"
              f"  ({mapping[cs]['name']})")
    if "oxygen" in mapping:
        print(f"  Oxygen exchange  : {mapping['oxygen']['exchange']}")
    if "co2" in mapping:
        print(f"  CO2 exchange     : {mapping['co2']['exchange']}")

    calc = TheoreticalYieldCalculator(
        metabolite_mapping=mapping,
        uptake_rate=10.0,           # mmol/gDW/h substrate uptake rate
        biomass_reaction=BIOMASS_ID,
    )

    # ── Context-specific model ──────────────────────────────────────────────
    try:
        ctx_model = results_obj.get_model(sample)
    except Exception as e:
        print(f"  [!] get_model failed ({e}), using bound-change model")
        rd        = results_obj.reaction_matrix.loc[sample].to_dict()
        ctx_model = recon._model.copy()
        for r in ctx_model.reactions:
            if not rd.get(r.id, False):
                r.bounds = (0.0, 0.0)
    ctx_model.objective = BIOMASS_ID

    # ── tINIT model yields ──────────────────────────────────────────────────
    print(f"\n  [Context-specific model — tINIT {sample}]")
    ctx_results = calc.calculate_yields(
        ctx_model,
        carbon_sources=carbon_sources,
        aerobic=True,
        anaerobic=True,
        product="biomass",
    )

    print(f"  {'Carbon Source':<14} {'Condition':<12} {'Yield (gDW/mmol)':>18}  {'Status':>6}")
    print(f"  {'─'*54}")
    for r in ctx_results:
        cond  = "Aerobic" if r.aerobic else "Anaerobic"
        y_str = f"{r.theoretical_yield:.4f}" if r.feasible else "N/A"
        ok    = "[OK]" if r.feasible else "[!] "
        print(f"  {ok} {r.carbon_source:<12} {cond:<12} {y_str:>18}")

    # ── Reference model yields (aerobic only) ───────────────────────────────
    print(f"\n  [Reference model — Human-GEM_consistent (aerobic only)]")
    ref_results = calc.calculate_yields(
        recon._model,
        carbon_sources=carbon_sources,
        aerobic=True,
        anaerobic=False,
        product="biomass",
    )
    print(f"  {'Carbon Source':<14} {'Condition':<12} {'Yield (gDW/mmol)':>18}  {'Status':>6}")
    print(f"  {'─'*54}")
    for r in ref_results:
        cond  = "Aerobic" if r.aerobic else "Anaerobic"
        y_str = f"{r.theoretical_yield:.4f}" if r.feasible else "N/A"
        ok    = "[OK]" if r.feasible else "[!] "
        print(f"  {ok} {r.carbon_source:<12} {cond:<12} {y_str:>18}")

    # ── Detailed comparison table ───────────────────────────────────────────
    ctx_df = calc.results_to_dataframe(ctx_results)
    ref_df = calc.results_to_dataframe(ref_results)
    ctx_df["model"] = f"tINIT_{sample}"
    ref_df["model"] = "reference"
    combined = pd.concat([ref_df, ctx_df], ignore_index=True)

    # Pivot: rows = carbon source + condition, cols = reference vs tINIT
    try:
        pivot = combined.pivot_table(
            index=["carbon_source", "aerobic"],
            columns="model",
            values="theoretical_yield",
            aggfunc="first",
        )
        print(f"\n  [Comparison: reference vs tINIT theoretical yields]")
        print(f"  {pivot.to_string()}")
    except Exception:
        pass

    # ── Save CSV ────────────────────────────────────────────────────────────
    os.makedirs(save_dir, exist_ok=True)
    out_csv = os.path.join(save_dir, "theoretical_yields.csv")
    combined.to_csv(out_csv, index=False)
    print(f"\n  [OK] yields saved -> {out_csv}")

    # ── Yield visualization ─────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from troppo.visualization import plot_yield_comparison

        fig = plot_yield_comparison(
            results_obj, ctx_df,
            save_path=os.path.join(save_dir, "theoretical_yields.png"),
        )
        plt.close(fig)
        print(f"  [OK] yield chart -> theoretical_yields.png")
    except Exception as e:
        print(f"  [!!] yield chart -> {e}")

    return combined


# ── Main ───────────────────────────────────────────────────────────────────────
def main(args):
    print(f"\n{SEP}")
    print("  MCF7 Context-Specific GEM — tINIT + Validation + Yields")
    print(SEP)
    print(f"  Model      : {MODEL_PATH}")
    print(f"  Expression : {EXPR_PATH}")
    print(f"  Tasks      : {TASK_PATH}")
    print(f"  Sample     : {args.sample}")
    print(f"  Solver     : {args.solver}")

    # ── Load expression data ───────────────────────────────────────────────────
    expr_df = load_ccle(EXPR_PATH)
    print(f"\n  Expression : {expr_df.shape[0]} genes x {expr_df.shape[1]} samples")
    if args.sample not in expr_df.columns:
        raise ValueError(
            f"Sample '{args.sample}' not found.\n"
            f"Available: {expr_df.columns.tolist()}"
        )

    # ── Initialize TINITReconstructor ──────────────────────────────────────────
    recon = TINITReconstructor(
        model_path=MODEL_PATH,
        expression_data=expr_df,
        gene_id_type="ensembl_gene_id",
        solver=args.solver,
        verbose=False,
    )
    print(f"  Model      : {len(recon._model.reactions)} rxns,"
          f" {len(recon._model.metabolites)} mets")

    tasks = recon._load_tasks(TASK_PATH)
    tasks = recon._filter_compatible_tasks(tasks)
    print(f"  Tasks      : {len(tasks)} compatible essential tasks")

    os.makedirs(RESULTS_DIR, exist_ok=True)

    # ── Stage 1: Reference evaluation ─────────────────────────────────────────
    ref_results = stage1_reference(recon, tasks)
    ref_passed  = sum(ref_results.values())
    if ref_passed < len(ref_results):
        print(f"\n  [!] Reference: {len(ref_results) - ref_passed} task(s) failed")
    else:
        print(f"\n  [OK] Reference: {ref_passed}/{len(ref_results)} tasks passed")

    # ── Stage 2: FVA essential reaction detection ──────────────────────────────
    essential, fva_df = stage2_fva_essential(recon)
    fva_df.to_csv(
        os.path.join(RESULTS_DIR, "fva_essential_reactions.csv"), index=False
    )

    # ── Stage 3: tINIT reconstruction ─────────────────────────────────────────
    results_obj = stage3_reconstruction(recon, args.sample, essential)

    # ── Save reconstructed model (SBML + reaction matrix CSV + pickle) ─────────
    print(f"\n{SEP}")
    print("  Saving Reconstructed Model")
    print(SEP)

    # SBML (.xml) — one file per sample, named <sample>.xml
    try:
        results_obj.to_sbml(RESULTS_DIR)
        xml_path = os.path.join(RESULTS_DIR, f"{args.sample}.xml")
        print(f"  [OK] SBML model        -> {args.sample}.xml")
    except Exception as e:
        print(f"  [!!] SBML save failed  -> {e}")

    # Reaction activity matrix (.csv) — rows=samples, cols=reactions, values=bool
    try:
        csv_path = os.path.join(RESULTS_DIR, "reaction_matrix.csv")
        results_obj.to_csv(csv_path)
        print(f"  [OK] Reaction matrix   -> reaction_matrix.csv")
    except Exception as e:
        print(f"  [!!] Reaction matrix   -> {e}")

    # Serialised results object (.pkl) — can be reloaded with ReconstructionResults.load()
    try:
        pkl_path = os.path.join(RESULTS_DIR, "reconstruction_results.pkl")
        results_obj.save(pkl_path)
        print(f"  [OK] Results pickle    -> reconstruction_results.pkl")
    except Exception as e:
        print(f"  [!!] Pickle save       -> {e}")

    # ── Stage 4: Visualization ─────────────────────────────────────────────────
    if not args.no_plots:
        stage4_visualize(results_obj, RESULTS_DIR)
    else:
        print(f"\n  [Stage 4 skipped — --no-plots]")

    # ── Stage 5: Biomass flux ──────────────────────────────────────────────────
    v_ref, v_tinit = stage5_biomass(recon, results_obj, args.sample)

    # ── Stage 6: Theoretical yields ────────────────────────────────────────────
    yield_df = stage6_yields(recon, results_obj, args.sample, RESULTS_DIR)

    # ── Final summary ──────────────────────────────────────────────────────────
    task_res = results_obj._task_results.get(args.sample)
    t_pass   = sum(1 for v in task_res.results.values() if v) if task_res else 0
    t_total  = len(task_res.results) if task_res else 0

    print(f"\n{SEP}")
    print("  FINAL SUMMARY")
    print(SEP)
    print(f"  {'Metric':<44} {'Value':>16}")
    print(f"  {'─'*62}")
    print(f"  {'Reference essential tasks passed':<44} {ref_passed:>3}/{len(ref_results):<12}")
    print(f"  {'tINIT essential tasks passed':<44} {t_pass:>3}/{t_total:<12}")
    print(f"  {'FVA essential reactions forced into tINIT':<44} {len(essential):>16}")
    print(f"  {'Reference biomass flux  (mmol/gDW/h)':<44} {v_ref:>16.4f}")
    print(f"  {'tINIT biomass flux  (mmol/gDW/h)':<44} {v_tinit:>16.4f}")
    if v_ref > 1e-9:
        print(f"  {'Biomass retention (%)':<44} {v_tinit/v_ref*100:>16.1f}")

    if not yield_df.empty:
        tinit_yields = yield_df[yield_df["model"].str.startswith("tINIT")]
        glu_aero = tinit_yields[
            (tinit_yields["carbon_source"] == "glucose") &
            (tinit_yields["aerobic"] == True)
        ]
        if not glu_aero.empty and glu_aero.iloc[0]["feasible"]:
            y = glu_aero.iloc[0]["theoretical_yield"]
            print(f"  {'Glucose aerobic yield  (gDW/mmol)':<44} {y:>16.4f}")

    print(f"  {'─'*62}")
    print(f"\n  Results directory: {RESULTS_DIR}")
    print(f"  Saved files:")
    print(f"    {args.sample}.xml             — context-specific model (SBML)")
    print(f"    reaction_matrix.csv           — boolean reaction activity matrix")
    print(f"    reconstruction_results.pkl    — serialised ReconstructionResults")
    print(f"    fva_essential_reactions.csv   — FVA result (all reactions)")
    print(f"    theoretical_yields.csv        — reference + tINIT yields")
    if not args.no_plots:
        print(f"    model_sizes.png               — active reactions per sample")
        print(f"    subsystem_coverage.png        — subsystem activation heatmap")
        print(f"    reaction_frequency.png        — reaction activation frequency")
        print(f"    theoretical_yields.png        — yield comparison bar chart")
    print(SEP)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "MCF7 context-specific GEM reconstruction via tINIT\n"
            "Includes essential task validation, visualization, "
            "biomass flux, and theoretical yield calculation."
        )
    )
    parser.add_argument(
        "--sample", default="global_0_nan_nan",
        help="CCLE sample name (default: global_0_nan_nan)",
    )
    parser.add_argument(
        "--solver", default="GUROBI", choices=["GUROBI", "CPLEX"],
        help="LP/MILP solver (default: GUROBI)",
    )
    parser.add_argument(
        "--no-plots", dest="no_plots", action="store_true",
        help="Skip visualization — useful in headless / no-display environments",
    )
    args = parser.parse_args()
    main(args)
