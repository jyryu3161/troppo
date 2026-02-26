"""
Essential Task Reaction-Level Analysis: Reference vs tINIT
==========================================================
Evaluates each essential task in the Human-GEM_consistent.xml reference model
and in a tINIT-reconstructed model, then reports per-category and per-reaction
pass/fail results in detail.

Design principles
-----------------
- Tasks with metabolites absent from the model are not removed; they are recorded
  as FAIL(incompatible) so that pass rates remain honest.
  → classify_tasks separates (classifies) tasks but includes all in the output.
- Reference and tINIT results are written side by side in the same Excel file.

Output file: results/essential_task_detail/task_reaction_detail.xlsx
  summary          : category-level pass rates (ref vs tINIT)
  task_level       : per-task pass/fail + status (OK/LOST/GAINED/REF_FAIL/INCOMPATIBLE)
  reaction_detail  : per-reaction/metabolite required bounds + ref/tINIT actual flux comparison

Usage
-----
  python analyze_biomass.py [--sample SAMPLE_NAME]
  e.g. python analyze_biomass.py --sample global_0_nan_nan
"""

import os, sys, re, warnings, argparse
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

import cobra
import pandas as pd
import numpy as np

from troppo.tasks.task_io import ExcelTaskIO
from troppo.tasks.core import TaskEvaluator

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH  = os.path.join(BASE_DIR, "..", "..", "tests", "data", "Human-GEM.xml")
TASK_PATH   = os.path.join(BASE_DIR, "..", "..", "tests", "data", "metabolicTasks_Essential.xlsx")
EXPR_PATH   = os.path.join(BASE_DIR, "data", "CCLE_expression_ACH-000019_scores.csv")
RESULTS_DIR = os.path.join(BASE_DIR, "results", "essential_task_detail")

CATEGORY_DESCRIPTIONS = {
    "ER": "Energy & Redox (ATP rephosphorylation, TCA, OXPHOS)",
    "BS": "Biosynthesis (nucleotides, cholesterol, phospholipids, vitamins)",
    "SU": "Substrate Utilization (amino acids, fatty acids, choline, inositol)",
    "IC": "Internal Conversions (glycolysis, acetyl-CoA, TCA intermediates)",
    "GR": "Growth (biomass production on Ham's F12 medium)",
}
_ALLMETSIN = {'ALLMETSIN[s]', 'ALLMETSIN[e]', 'ALLMETSIN'}


# ── Biomass flux calculation ───────────────────────────────────────────────────

def get_biomass_reaction_id(model: cobra.Model) -> str:
    """Return the objective (biomass) reaction ID from the model."""
    try:
        from cobra.util.solver import linear_reaction_coefficients
        coeffs = linear_reaction_coefficients(model)
        if coeffs:
            return max(coeffs, key=lambda r: abs(coeffs[r])).id
    except Exception:
        pass
    for rxn in model.reactions:
        if 'biomass' in rxn.id.lower():
            return rxn.id
    return None


def calculate_biomass_flux(model: cobra.Model,
                           biomass_rxn_id: str,
                           bound_changes: dict = None) -> dict:
    """
    Calculate biomass flux via FBA.

    Parameters
    ----------
    model           : COBRApy model (original bounds protected via context)
    biomass_rxn_id  : biomass reaction ID
    bound_changes   : {rxn_id: (lb, ub)} — inactive reaction bounds from tINIT

    Returns
    -------
    dict: status, objective_value, biomass_flux
    """
    with model:
        if bound_changes:
            for rxn_id, (lb, ub) in bound_changes.items():
                try:
                    model.reactions.get_by_id(rxn_id).bounds = (lb, ub)
                except Exception:
                    pass
        sol = model.optimize()
        if sol.status == 'optimal':
            obj_val  = sol.objective_value
            bio_flux = sol.fluxes.get(biomass_rxn_id, float('nan'))
        else:
            obj_val  = float('nan')
            bio_flux = float('nan')
    return {
        "status":          sol.status,
        "objective_value": obj_val,
        "biomass_flux":    bio_flux,
    }


# ── FVA-based biomass essential reaction extraction ───────────────────────────

def get_biomass_essential_reactions(model: cobra.Model,
                                    biomass_rxn_id: str,
                                    min_fraction: float = 0.01,
                                    tol: float = 1e-6) -> tuple:
    """
    Extract biomass-essential reactions via FVA.

    Runs FVA on the reference model while constraining
    biomass flux >= min_fraction * optimal, and returns reactions
    that must carry non-zero flux.

    Parameters
    ----------
    model           : COBRApy model (original bounds protected via context)
    biomass_rxn_id  : biomass reaction ID
    min_fraction    : minimum fraction of optimal biomass to maintain (default 1%)
    tol             : tolerance for considering flux non-zero

    Returns
    -------
    (essential_rxn_ids : list[str], fva_df : pd.DataFrame)
        fva_df columns: rxn_id, min_flux, max_flux, is_essential, rxn_name, rxn_equation
    """
    with model:
        # 1. Optimize biomass
        model.objective = biomass_rxn_id
        opt_sol = model.optimize()

        if opt_sol.status != 'optimal' or opt_sol.objective_value < tol:
            print(f"  [WARNING] Reference model cannot produce biomass "
                  f"(status: {opt_sol.status})")
            return [], pd.DataFrame()

        opt_biomass = opt_sol.objective_value
        min_biomass = opt_biomass * min_fraction
        print(f"  Optimal biomass  : {opt_biomass:.6f}")
        print(f"  Min constraint   : {min_biomass:.8f} "
              f"({min_fraction * 100:.1f}% of optimal)")

        # 2. Enforce minimum biomass lower bound
        model.reactions.get_by_id(biomass_rxn_id).lower_bound = min_biomass

        # 3. Run FVA
        #    fraction_of_optimum=0.0: the biomass lb already constrains feasibility,
        #    so the full flux range is explored without an additional objective constraint.
        print(f"  Running FVA on {len(model.reactions)} reactions "
              f"(processes=1, may take a few minutes) ...")
        fva = cobra.flux_analysis.flux_variability_analysis(
            model,
            reaction_list=None,
            fraction_of_optimum=0.0,
            processes=1,
        )

    # 4. Extract essential reactions
    #    min_flux > tol  : must carry forward flux
    #    max_flux < -tol : must carry reverse flux
    essential = []
    rows = []
    for rxn_id, row in fva.iterrows():
        min_f = row['minimum']
        max_f = row['maximum']
        is_ess = (min_f > tol) or (max_f < -tol)
        if is_ess:
            essential.append(rxn_id)
        try:
            rxn      = model.reactions.get_by_id(rxn_id)
            rxn_name = rxn.name
            rxn_eq   = rxn.reaction
        except Exception:
            rxn_name = ''
            rxn_eq   = ''
        rows.append({
            'rxn_id':       rxn_id,
            'min_flux':     round(min_f, 8),
            'max_flux':     round(max_f, 8),
            'is_essential': is_ess,
            'rxn_name':     rxn_name,
            'rxn_equation': rxn_eq,
        })

    fva_df = pd.DataFrame(rows)
    return essential, fva_df


# ── Utilities ──────────────────────────────────────────────────────────────────

def get_category(task_name: str) -> str:
    for cat in CATEGORY_DESCRIPTIONS:
        if task_name.startswith(cat):
            return cat
    return "OTHER"


def translate_task_metabolite_ids(tasks: list, model: cobra.Model) -> list:
    """Translate task metabolite IDs (name[compartment]) to model HMDB IDs."""
    name_comp_to_id = {}
    for m in model.metabolites:
        name_comp_to_id[f"{m.name.lower()}[{m.compartment}]"] = m.id

    def tr(k):
        if k in _ALLMETSIN:
            return k
        return name_comp_to_id.get(k.lower(), k)

    for task in tasks:
        task.inflow_dict      = {tr(k): v for k, v in task.inflow_dict.items()}
        task.outflow_dict     = {tr(k): v for k, v in task.outflow_dict.items()}
        task.flux_constraints = {tr(k): v for k, v in task.flux_constraints.items()}
        task.reaction_dict    = {rk: ({tr(m): c for m, c in st.items()}, bnd)
                                 for rk, (st, bnd) in task.reaction_dict.items()}
    return tasks


def classify_tasks(tasks: list, model: cobra.Model):
    """
    Classify tasks as compatible / incompatible.

    Incompatible tasks are not removed; they are returned as
    {task_name: [missing_mets]} so they can be recorded as FAIL(incompatible).
    """
    model_met_ids = set(m.id for m in model.metabolites)
    compatible, incompatible = [], {}
    for task in tasks:
        real_mets = [m for m in list(task.inflow_dict) + list(task.outflow_dict)
                     if m not in _ALLMETSIN]
        missing = [m for m in real_mets if m not in model_met_ids]
        if missing:
            incompatible[task.name] = missing
        else:
            compatible.append(task)
    return compatible, incompatible


def load_ccle_expression(path: str) -> pd.DataFrame:
    """Load CCLE scores CSV as a genes × samples DataFrame."""
    df = pd.read_csv(path)
    gene_cols = [c for c in df.columns if "ENSG" in str(c)]
    meta_cols = [c for c in df.columns if c.startswith("Unnamed")]
    df_expr   = df[gene_cols].copy()
    df_expr.index = df[meta_cols].apply(
        lambda r: "_".join(str(v) for v in r), axis=1)
    df_T = df_expr.T
    df_T.index = [re.search(r"(ENSG\d+)", str(c)).group(1)
                  if re.search(r"(ENSG\d+)", str(c)) else c
                  for c in df_T.index]
    return df_T


def get_reaction_equation(model: cobra.Model, rxn_id: str) -> str:
    try:
        return model.reactions.get_by_id(rxn_id).reaction
    except Exception:
        return ""


def get_met_name(model: cobra.Model, met_id: str) -> str:
    try:
        return model.metabolites.get_by_id(met_id).name
    except Exception:
        return met_id


# ── Core: collect detail rows ─────────────────────────────────────────────────

def collect_detail_rows(tasks: list,
                        incompatible: dict,
                        raw_ref: dict,
                        raw_tinit: dict,
                        cobamp_rxn_names: list,
                        model: cobra.Model) -> list:
    """
    Collect reaction_detail rows for all tasks (compatible + incompatible).

    Parameters
    ----------
    tasks        : compatible task list (evaluated by TaskEvaluator)
    incompatible : {task_name: [missing_met_ids]}
    raw_ref      : batch_evaluate result [0] — reference
    raw_tinit    : batch_evaluate result [1] — tINIT
    """

    def get_flux(sol, rxn_name: str) -> float:
        if sol is None:
            return float('nan')
        try:
            return float(sol[rxn_name])
        except Exception:
            pass
        try:
            return float(sol[cobamp_rxn_names.index(rxn_name)])
        except Exception:
            return float('nan')

    def flux_ok_str(flux, lb, ub=None, is_outflow=False) -> str:
        if np.isnan(flux):
            return "N/A"
        if is_outflow:
            return "YES" if flux >= lb - 1e-6 else "NO"
        return "YES" if lb - 1e-6 <= flux <= (ub if ub is not None else 1e18) + 1e-6 else "NO"

    rows = []

    # ── Compatible tasks ───────────────────────────────────────────────────────
    for task in tasks:
        ref_truth,   _, sol_ref   = raw_ref.get(task.name,   (False, {}, None))
        tinit_truth, _, sol_tinit = raw_tinit.get(task.name, (False, {}, None))
        cat        = get_category(task.name)
        ref_str    = "PASS" if ref_truth   else "FAIL"
        tinit_str  = "PASS" if tinit_truth else "FAIL"

        def make_row(item_type, item_id, item_name, rxn_id, rxn_eq,
                     lb, ub, note="", is_outflow=False):
            f_ref   = get_flux(sol_ref,   rxn_id)
            f_tinit = get_flux(sol_tinit, rxn_id)
            return {
                "task_id":       task.name,
                "category":      cat,
                "ref_result":    ref_str,
                "tinit_result":  tinit_str,
                "item_type":     item_type,
                "item_id":       item_id,
                "item_name":     item_name,
                "reaction_id":   rxn_id,
                "reaction_eq":   rxn_eq,
                "required_lb":   lb,
                "required_ub":   ub,
                "ref_flux":      round(f_ref,   6) if not np.isnan(f_ref)   else "",
                "ref_flux_ok":   flux_ok_str(f_ref,   lb, ub, is_outflow),
                "tinit_flux":    round(f_tinit, 6) if not np.isnan(f_tinit) else "",
                "tinit_flux_ok": flux_ok_str(f_tinit, lb, ub, is_outflow),
                "note":          note,
            }

        for met_id, (lb, ub) in task.inflow_dict.items():
            if met_id in _ALLMETSIN:
                continue
            rows.append(make_row(
                "INFLOW", met_id, get_met_name(model, met_id),
                f"{met_id}_inflow", "", lb, ub))

        for met_id, (lb, ub) in task.outflow_dict.items():
            if met_id in _ALLMETSIN:
                continue
            note = "biomass_output" if ("biomass" in met_id.lower()
                                        or "temp001" in met_id.lower()) else ""
            rows.append(make_row(
                "OUTFLOW", met_id, get_met_name(model, met_id),
                f"{met_id}_outflow", "", lb, ub, note, is_outflow=True))

        for rxn_id, (lb, ub) in task.flux_constraints.items():
            rows.append(make_row(
                "FLUX_CONSTRAINT", rxn_id, "",
                rxn_id, get_reaction_equation(model, rxn_id),
                lb, ub, "existing_rxn_bounds_changed"))

        for rxn_key, (stoich, (lb, ub)) in task.reaction_dict.items():
            rxn_name = f"{task.name}_{rxn_key}_task_reaction"
            reactants = [f"{abs(v)} {k}" for k, v in stoich.items() if v < 0]
            products  = [f"{abs(v)} {k}" for k, v in stoich.items() if v > 0]
            eq_str    = " + ".join(reactants) + " => " + " + ".join(products)
            rows.append(make_row(
                "CUSTOM_REACTION", rxn_key, "",
                rxn_name, eq_str, lb, ub, "task-specific reaction"))

    # ── Incompatible tasks → record as FAIL(incompatible) ────────────────────
    for task_name, missing_mets in incompatible.items():
        cat = get_category(task_name)
        rows.append({
            "task_id":       task_name,
            "category":      cat,
            "ref_result":    "FAIL",
            "tinit_result":  "FAIL",
            "item_type":     "INCOMPATIBLE",
            "item_id":       "; ".join(missing_mets),
            "item_name":     "missing from model",
            "reaction_id":   "",
            "reaction_eq":   "",
            "required_lb":   "",
            "required_ub":   "",
            "ref_flux":      "",
            "ref_flux_ok":   "N/A",
            "tinit_flux":    "",
            "tinit_flux_ok": "N/A",
            "note":          "metabolite_not_in_model → counted as FAIL",
        })

    return rows


# ── Result collection ─────────────────────────────────────────────────────────

def collect_task_level(all_tasks_ordered: list,
                       compatible: list,
                       incompatible: dict,
                       raw_ref: dict,
                       raw_tinit: dict) -> list:
    """Collect task_level sheet rows."""
    comp_names = {t.name for t in compatible}
    rows = []
    for task in all_tasks_ordered:
        cat = get_category(task.name)
        ann = task.annotations if hasattr(task, "annotations") and task.annotations else {}
        desc = ann.get("DESCRIPTION", ann.get("description", ""))
        is_incompatible = task.name in incompatible

        if is_incompatible:
            ref_truth   = False
            tinit_truth = False
            status      = "INCOMPATIBLE"
        else:
            ref_truth,   _, _ = raw_ref.get(task.name,   (False, {}, None))
            tinit_truth, _, _ = raw_tinit.get(task.name, (False, {}, None))
            if ref_truth and tinit_truth:
                status = "OK"
            elif ref_truth and not tinit_truth:
                status = "LOST"
            elif not ref_truth and tinit_truth:
                status = "GAINED"
            else:
                status = "REF_FAIL"

        rows.append({
            "task_id":            task.name,
            "category":           cat,
            "category_desc":      CATEGORY_DESCRIPTIONS.get(cat, ""),
            "description":        desc,
            "should_fail":        task.should_fail,
            "incompatible":       is_incompatible,
            "incompatible_reason": ("; ".join(incompatible[task.name])
                                    if is_incompatible else ""),
            "ref_result":         "FAIL(incompatible)" if is_incompatible else
                                  ("PASS" if ref_truth else "FAIL"),
            "ref_pass":           ref_truth,
            "tinit_result":       "FAIL(incompatible)" if is_incompatible else
                                  ("PASS" if tinit_truth else "FAIL"),
            "tinit_pass":         tinit_truth,
            "status":             status,
            "n_inflow":           len([k for k in task.inflow_dict if k not in _ALLMETSIN]),
            "n_outflow":          len([k for k in task.outflow_dict if k not in _ALLMETSIN]),
            "n_flux_constraints": len(task.flux_constraints),
            "n_custom_reactions": len(task.reaction_dict),
        })
    return rows


def collect_summary(task_rows: list) -> list:
    """Collect category-level pass rate summary rows (ref + tINIT)."""
    cat_stats = {cat: {"ref_pass":0,"tinit_pass":0,"total":0}
                 for cat in list(CATEGORY_DESCRIPTIONS)+["OTHER"]}
    for r in task_rows:
        cat = r["category"]
        cat_stats[cat]["total"]      += 1
        cat_stats[cat]["ref_pass"]   += int(r["ref_pass"])
        cat_stats[cat]["tinit_pass"] += int(r["tinit_pass"])

    rows = []
    for cat, desc in CATEGORY_DESCRIPTIONS.items():
        s = cat_stats[cat]
        if s["total"] == 0:
            continue
        rows.append({
            "category":       cat,
            "description":    desc,
            "ref_pass":       s["ref_pass"],
            "ref_total":      s["total"],
            "ref_rate_%":     round(s["ref_pass"]   / s["total"] * 100, 1),
            "ref_fail":       s["total"] - s["ref_pass"],
            "tinit_pass":     s["tinit_pass"],
            "tinit_total":    s["total"],
            "tinit_rate_%":   round(s["tinit_pass"] / s["total"] * 100, 1),
            "tinit_fail":     s["total"] - s["tinit_pass"],
        })
    total = sum(s["total"]      for s in cat_stats.values())
    rp    = sum(s["ref_pass"]   for s in cat_stats.values())
    tp    = sum(s["tinit_pass"] for s in cat_stats.values())
    rows.append({
        "category":    "TOTAL",
        "description": "All essential tasks",
        "ref_pass":    rp,   "ref_total":   total, "ref_rate_%":   round(rp/max(total,1)*100,1),
        "ref_fail":    total-rp,
        "tinit_pass":  tp,   "tinit_total": total, "tinit_rate_%": round(tp/max(total,1)*100,1),
        "tinit_fail":  total-tp,
    })
    return rows


# ── Excel output ───────────────────────────────────────────────────────────────

def write_excel(out_path: str, sheets: dict):
    """sheets = {sheet_name: DataFrame}"""
    with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
        for sname, df in sheets.items():
            df.to_excel(writer, sheet_name=sname, index=False)
            ws = writer.sheets[sname]
            for col_idx, col in enumerate(df.columns, 1):
                max_len = max(len(str(col)),
                              df[col].astype(str).str.len().max() if len(df) else 0)
                ws.column_dimensions[
                    ws.cell(1, col_idx).column_letter
                ].width = min(max_len + 2, 60)


# ── Main ───────────────────────────────────────────────────────────────────────

def main(sample_name: str):
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # 1. Load model
    print("=== 1. Loading Human-GEM.xml (tests/data) ===")
    model = cobra.io.read_sbml_model(MODEL_PATH)
    print(f"  Reactions: {len(model.reactions)}, Metabolites: {len(model.metabolites)}")

    # 2. Load and translate tasks
    print("\n=== 2. Loading & translating tasks ===")
    tasks_all = ExcelTaskIO().read_task(TASK_PATH)
    tasks_all = translate_task_metabolite_ids(tasks_all, model)
    compatible, incompatible = classify_tasks(tasks_all, model)
    print(f"  Total tasks   : {len(tasks_all)}")
    print(f"  Compatible    : {len(compatible)}")
    print(f"  Incompatible  : {len(incompatible)}  ← recorded as FAIL instead of removed")
    for tid, miss in incompatible.items():
        print(f"    [INCOMPATIBLE] {tid}: {miss}")

    # 3a+3b. tINIT reconstruction (FVA essential extraction → reconstruction)
    print(f"\n=== 3a. TINITReconstructor initialization ===")
    from troppo import TINITReconstructor
    expr_df = load_ccle_expression(EXPR_PATH)
    if sample_name not in expr_df.columns:
        avail = expr_df.columns.tolist()
        raise ValueError(f"Sample '{sample_name}' not found. Available: {avail}")
    print(f"  Expression: {expr_df.shape[0]} genes × {expr_df.shape[1]} samples")

    recon = TINITReconstructor(
        model_path=os.path.abspath(MODEL_PATH),
        expression_data=expr_df,
        gene_id_type="ensembl_gene_id",
        solver="GUROBI",
        verbose=False,
    )

    print(f"\n=== 3b. FVA-based Biomass Essential Reaction Extraction ===")
    biomass_rxn_id = get_biomass_reaction_id(model)
    fva_essential, fva_df = recon.derive_biomass_essential_reactions(
        biomass_rxn_id=biomass_rxn_id,
        min_fraction=0.01,
    )

    print(f"\n=== 3c. tINIT reconstruction (sample: {sample_name}) ===")
    recon.configure(essential_reactions=fva_essential)
    print(f"  essential_reactions (FVA-derived): {len(fva_essential)} reactions")
    results_obj = recon.run(samples=[sample_name], gap_fill=False, validate_tasks=False)
    rxn_matrix  = results_obj.reaction_matrix
    result_dict = rxn_matrix.loc[sample_name].to_dict()
    n_active    = sum(1 for v in result_dict.values() if v)
    print(f"  Active reactions: {n_active}/{len(model.reactions)}")

    # tINIT bound_changes: inactive reactions → (0, 0)
    tinit_bc = {rxn.id: (0.0, 0.0)
                for rxn in model.reactions
                if not result_dict.get(rxn.id, False)}
    print(f"  Inactivated reactions: {len(tinit_bc)}")

    # 3d. Biomass flux (FBA)
    print("\n=== 3d. Biomass Flux Calculation (FBA) ===")
    print(f"  Biomass reaction: {biomass_rxn_id}")

    ref_bio   = calculate_biomass_flux(model, biomass_rxn_id, {})
    tinit_bio = calculate_biomass_flux(model, biomass_rxn_id, tinit_bc)

    print(f"  Reference → status: {ref_bio['status']},  "
          f"biomass flux: {ref_bio['biomass_flux']:.6f}")
    print(f"  tINIT     → status: {tinit_bio['status']},  "
          f"biomass flux: {tinit_bio['biomass_flux']:.6f}")
    if (not np.isnan(ref_bio['biomass_flux'])
            and ref_bio['biomass_flux'] > 1e-9
            and not np.isnan(tinit_bio['biomass_flux'])):
        ratio = tinit_bio['biomass_flux'] / ref_bio['biomass_flux'] * 100
        print(f"  tINIT / Reference: {ratio:.1f}%")

    # 4. TaskEvaluator — evaluate reference ({}) and tINIT (tinit_bc) simultaneously
    print("\n=== 4. TaskEvaluator (ref + tINIT, output_sol=True) ===")
    tev = TaskEvaluator(model=model, tasks=compatible, solver='GUROBI')
    raw = tev.batch_evaluate(
        bound_changes=[{}, tinit_bc],  # index 0 = ref, index 1 = tINIT
        threads=1,
        output_sol=True,
    )
    raw_ref, raw_tinit = raw[0], raw[1]
    cobamp_rxn_names   = list(tev.model.reaction_names)
    print("  Done.")

    # 5. Collect results
    print("\n=== 5. Collecting results ===")
    task_rows   = collect_task_level(tasks_all, compatible, incompatible,
                                     raw_ref, raw_tinit)
    summary_rows = collect_summary(task_rows)
    detail_rows  = collect_detail_rows(compatible, incompatible,
                                       raw_ref, raw_tinit,
                                       cobamp_rxn_names, model)

    summary_df = pd.DataFrame(summary_rows)
    task_df    = pd.DataFrame(task_rows)
    detail_df  = pd.DataFrame(detail_rows)

    # 6. Excel output
    biomass_rows = [
        {
            "model":            "Reference (full model)",
            "sample":           "—",
            "biomass_reaction": biomass_rxn_id,
            "status":           ref_bio["status"],
            "biomass_flux":     ref_bio["biomass_flux"],
            "ratio_%":          100.0,
        },
        {
            "model":            f"tINIT",
            "sample":           sample_name,
            "biomass_reaction": biomass_rxn_id,
            "status":           tinit_bio["status"],
            "biomass_flux":     tinit_bio["biomass_flux"],
            "ratio_%":          (
                round(tinit_bio["biomass_flux"] / ref_bio["biomass_flux"] * 100, 2)
                if ref_bio["biomass_flux"] > 1e-9
                   and not np.isnan(tinit_bio["biomass_flux"])
                else float('nan')
            ),
        },
    ]
    biomass_df = pd.DataFrame(biomass_rows)

    out_path = os.path.join(RESULTS_DIR, "task_reaction_detail.xlsx")
    write_excel(out_path, {
        "biomass_flux":    biomass_df,
        "fva_essential":   fva_df,
        "summary":         summary_df,
        "task_level":      task_df,
        "reaction_detail": detail_df,
    })
    print(f"  Written: {out_path}")

    # 7. Console summary
    print("\n" + "=" * 72)
    print(f"  Essential Task Summary  |  Reference vs tINIT ({sample_name})")
    print("=" * 72)
    print(f"  {'Cat':<6}  {'Ref':>10}  {'tINIT':>10}  Description")
    print("  " + "-" * 70)
    for row in summary_rows:
        rpt = f"{row['ref_pass']}/{row['ref_total']} ({row['ref_rate_%']}%)"
        tpt = f"{row['tinit_pass']}/{row['tinit_total']} ({row['tinit_rate_%']}%)"
        print(f"  {row['category']:<6}  {rpt:>10}  {tpt:>10}  {row['description']}")

    for status, label in [
        ("LOST",          "▼ Ref PASS → tINIT FAIL (function lost)"),
        ("GAINED",        "▲ Ref FAIL → tINIT PASS (unexpected pass)"),
        ("REF_FAIL",      "✗ Ref FAIL (model limitation)"),
        ("INCOMPATIBLE",  "⊘ Metabolite absent → both FAIL"),
    ]:
        subset = task_df[task_df["status"] == status]
        if subset.empty:
            continue
        print(f"\n  {label} ({len(subset)}):")
        for _, r in subset.iterrows():
            reason = f"  [{r['incompatible_reason']}]" if r["incompatible"] else ""
            print(f"    {r['task_id']:<10} [{r['category']}] {r['description'][:45]}{reason}")

    print(f"\n  [Biomass Flux]  reaction: {biomass_rxn_id}")
    print(f"    Reference : {ref_bio['biomass_flux']:.6f}  (status: {ref_bio['status']})")
    print(f"    tINIT     : {tinit_bio['biomass_flux']:.6f}  (status: {tinit_bio['status']})")
    if ref_bio['biomass_flux'] > 1e-9 and not np.isnan(tinit_bio['biomass_flux']):
        print(f"    Ratio     : {tinit_bio['biomass_flux']/ref_bio['biomass_flux']*100:.1f}%")

    print(f"\n  Output: {out_path}")
    print("  Sheets: summary | task_level | reaction_detail")
    print(f"  reaction_detail rows: {len(detail_df)}")
    print(f"    INFLOW:           {(detail_df['item_type']=='INFLOW').sum()}")
    print(f"    OUTFLOW:          {(detail_df['item_type']=='OUTFLOW').sum()}")
    print(f"    FLUX_CONSTRAINT:  {(detail_df['item_type']=='FLUX_CONSTRAINT').sum()}")
    print(f"    CUSTOM_REACTION:  {(detail_df['item_type']=='CUSTOM_REACTION').sum()}")
    print(f"    INCOMPATIBLE:     {(detail_df['item_type']=='INCOMPATIBLE').sum()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", default="global_0_nan_nan",
                        help="CCLE sample name to use for tINIT reconstruction "
                             "(default: global_0_nan_nan)")
    args = parser.parse_args()
    main(args.sample)
