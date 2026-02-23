#!/usr/bin/env python
"""
Diagnose biomass=0 models and perform systematic gap-filling.

Strategy:
1. Load working model (Lung_Healthy, biomass=11.05) as reference
2. For each broken model, find biomass precursor dead-ends
3. Use iterative gap-filling from universal model to restore biomass
"""
import sys, os, time
import cobra
from cobra.io import read_sbml_model, write_sbml_model
import pandas as pd
import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODELS_DIR = os.path.join(BASE_DIR, "covid19_tinit_results", "models")
RESULTS_DIR = os.path.join(BASE_DIR, "covid19_tinit_results")

print("=" * 70)
print("Biomass=0 Diagnosis & Systematic Gap-Fill")
print("=" * 70)

# ── 1. Load universal model ──
print("\n[1] Loading universal model...")
t0 = time.time()
universal = read_sbml_model(os.path.join(BASE_DIR, "HumanGEM_Consistent_COVID19_HAM.xml"))
print(f"    Universal model: {len(universal.reactions)} reactions, {len(universal.metabolites)} metabolites")

# ── 2. Load reference working model (Lung_Healthy) ──
print("\n[2] Loading reference model (Lung_Healthy)...")
ref_model = read_sbml_model(os.path.join(MODELS_DIR, "Lung_Healthy_model.xml"))
ref_model.objective = 'biomass_human'
ref_sol = ref_model.optimize()
print(f"    Lung_Healthy: {len(ref_model.reactions)} rxns, biomass={ref_sol.objective_value:.6f}")
ref_rxn_ids = set(r.id for r in ref_model.reactions)

# ── 3. Identify biomass reaction and its precursors ──
print("\n[3] Analyzing biomass reaction precursors...")
biomass_rxn = universal.reactions.get_by_id('biomass_human')
# Reactants of biomass (consumed metabolites)
biomass_reactants = {m.id: m for m, coeff in biomass_rxn.metabolites.items() if coeff < 0}
biomass_products = {m.id: m for m, coeff in biomass_rxn.metabolites.items() if coeff > 0}
print(f"    Biomass reactants (precursors): {len(biomass_reactants)}")
print(f"    Biomass products: {len(biomass_products)}")

# ── 4. Analyze each broken model ──
broken_models = [
    "Heart_Healthy", "Heart_COVID19",
    "Liver_COVID19",
    "Kidney_Healthy",
    "Intestine_Healthy", "Intestine_COVID19"
]

print("\n[4] Diagnosing broken models...")
print("-" * 70)

all_diagnoses = []

for sample_name in broken_models:
    xml_path = os.path.join(MODELS_DIR, f"{sample_name}_model.xml")
    if not os.path.exists(xml_path):
        print(f"\n  {sample_name}: model file not found, skipping")
        continue

    sm = read_sbml_model(xml_path)
    if len(sm.reactions) == 0:
        print(f"\n  {sample_name}: empty model, skipping")
        continue

    sm.objective = 'biomass_human'
    sol = sm.optimize()
    bm = sol.objective_value if sol.objective_value is not None else 0.0
    sm_rxn_ids = set(r.id for r in sm.reactions)
    sm_met_ids = set(m.id for m in sm.metabolites)

    # Reactions in reference but not in this model
    missing_from_ref = ref_rxn_ids - sm_rxn_ids
    extra_vs_ref = sm_rxn_ids - ref_rxn_ids

    print(f"\n  {sample_name}: {len(sm.reactions)} rxns, biomass={bm:.6f}")
    print(f"    Missing vs Lung_Healthy: {len(missing_from_ref)} rxns")
    print(f"    Extra vs Lung_Healthy: {len(extra_vs_ref)} rxns")

    # Check each biomass precursor: can it be produced?
    blocked_precursors = []
    for met_id, met in biomass_reactants.items():
        if met_id not in sm_met_ids:
            blocked_precursors.append((met_id, met.name, "MISSING_METABOLITE"))
            continue
        # Check if any reaction can produce this metabolite
        sm_met = sm.metabolites.get_by_id(met_id)
        producers = [r for r in sm_met.reactions
                     if r.metabolites[sm_met] > 0 and r.id != 'biomass_human']
        consumers = [r for r in sm_met.reactions
                     if r.metabolites[sm_met] < 0 and r.id != 'biomass_human']
        # Also check reversible reactions that consume it (they can run backward)
        rev_producers = [r for r in consumers if r.reversibility]
        total_producers = len(producers) + len(rev_producers)
        if total_producers == 0:
            blocked_precursors.append((met_id, met.name, "NO_PRODUCER"))

    if blocked_precursors:
        print(f"    Blocked biomass precursors ({len(blocked_precursors)}):")
        for mid, mname, reason in blocked_precursors[:10]:
            print(f"      {mid}: {mname} [{reason}]")
        if len(blocked_precursors) > 10:
            print(f"      ... and {len(blocked_precursors)-10} more")

    all_diagnoses.append({
        'sample': sample_name,
        'n_reactions': len(sm.reactions),
        'biomass': bm,
        'missing_vs_ref': len(missing_from_ref),
        'blocked_precursors': len(blocked_precursors),
        'blocked_details': "; ".join(f"{mid}({reason})" for mid, _, reason in blocked_precursors)
    })

# ── 5. Systematic gap-fill using FBA-based approach ──
print("\n" + "=" * 70)
print("[5] Systematic Gap-Fill")
print("=" * 70)

gapfill_report = []

for sample_name in broken_models:
    xml_path = os.path.join(MODELS_DIR, f"{sample_name}_model.xml")
    if not os.path.exists(xml_path):
        continue

    sm = read_sbml_model(xml_path)
    if len(sm.reactions) == 0:
        print(f"\n  {sample_name}: empty model, skipping")
        continue

    sm.objective = 'biomass_human'
    sol = sm.optimize()
    bm = sol.objective_value if sol.objective_value is not None else 0.0

    if bm > 1e-6:
        print(f"\n  {sample_name}: biomass={bm:.6f} (already working)")
        continue

    print(f"\n  {sample_name}: biomass={bm:.6f}, attempting gap-fill...")

    sm_rxn_ids = set(r.id for r in sm.reactions)
    sm_met_ids = set(m.id for m in sm.metabolites)

    # Strategy: Add all reactions from reference model that are missing,
    # then prune back to minimal set needed for biomass.
    # But first, let's try a targeted approach: add reactions subsystem by subsystem.

    # Get reactions present in reference but not in this model
    candidate_rxns = []
    for rid in ref_rxn_ids:
        if rid not in sm_rxn_ids and rid in set(r.id for r in universal.reactions):
            candidate_rxns.append(rid)

    print(f"    Candidate gap-fill reactions (from Lung_Healthy): {len(candidate_rxns)}")

    # Group by subsystem for targeted analysis
    subsystem_rxns = {}
    for rid in candidate_rxns:
        urxn = universal.reactions.get_by_id(rid)
        ss = urxn.subsystem if urxn.subsystem else "Unknown"
        if ss not in subsystem_rxns:
            subsystem_rxns[ss] = []
        subsystem_rxns[ss].append(rid)

    # Try adding all candidates at once first to see if it works
    test_model = sm.copy()
    added_all = []
    test_met_ids = set(m.id for m in test_model.metabolites)
    test_rxn_ids = set(r.id for r in test_model.reactions)
    for rid in candidate_rxns:
        if rid not in test_rxn_ids:
            urxn = universal.reactions.get_by_id(rid).copy()
            for met in urxn.metabolites:
                if met.id not in test_met_ids:
                    test_model.add_metabolites([met.copy()])
                    test_met_ids.add(met.id)
            test_model.add_reactions([urxn])
            test_rxn_ids.add(rid)
            added_all.append(rid)

    test_model.objective = 'biomass_human'
    test_sol = test_model.optimize()
    test_bm = test_sol.objective_value if test_sol.objective_value is not None else 0.0
    print(f"    Adding ALL {len(added_all)} ref reactions → biomass={test_bm:.6f}")

    if test_bm < 1e-6:
        # Even with all reference reactions, no biomass. Try universal model reactions.
        print(f"    Reference reactions insufficient. Trying ALL universal model reactions...")
        test_model2 = sm.copy()
        added_all2 = []
        test_met_ids2 = set(m.id for m in test_model2.metabolites)
        test_rxn_ids2 = set(r.id for r in test_model2.reactions)
        for urxn in universal.reactions:
            if urxn.id not in test_rxn_ids2:
                urxn_copy = urxn.copy()
                for met in urxn_copy.metabolites:
                    if met.id not in test_met_ids2:
                        test_model2.add_metabolites([met.copy()])
                        test_met_ids2.add(met.id)
                test_model2.add_reactions([urxn_copy])
                test_rxn_ids2.add(urxn.id)
                added_all2.append(urxn.id)

        test_model2.objective = 'biomass_human'
        test_sol2 = test_model2.optimize()
        test_bm2 = test_sol2.objective_value if test_sol2.objective_value is not None else 0.0
        print(f"    Adding ALL {len(added_all2)} universal reactions → biomass={test_bm2:.6f}")

        if test_bm2 < 1e-6:
            print(f"    *** CRITICAL: Even universal model can't restore biomass!")
            print(f"    *** This suggests a fundamental stoichiometric issue.")
            # Check if biomass works in universal model itself
            universal.objective = 'biomass_human'
            u_sol = universal.optimize()
            print(f"    *** Universal model biomass: {u_sol.objective_value:.6f}")
            gapfill_report.append({
                'sample': sample_name,
                'status': 'UNFIXABLE',
                'reactions_added': 0,
                'biomass_after': 0.0,
                'details': 'Even full universal model cannot restore biomass'
            })
            continue
        else:
            # Universal works but reference doesn't - find minimal set from universal
            candidate_rxns = added_all2
            test_model = test_model2
            test_bm = test_bm2

    # Binary search for minimal gap-fill set
    print(f"    Finding minimal gap-fill set...")

    # First, identify which subsystems are needed
    needed_subsystems = set()
    for ss, rids in sorted(subsystem_rxns.items(), key=lambda x: len(x[1]), reverse=True):
        # Test removing this subsystem
        test_without = sm.copy()
        t_met_ids = set(m.id for m in test_without.metabolites)
        t_rxn_ids = set(r.id for r in test_without.reactions)
        for rid in candidate_rxns:
            urxn = universal.reactions.get_by_id(rid)
            curr_ss = urxn.subsystem if urxn.subsystem else "Unknown"
            if curr_ss == ss:
                continue  # skip this subsystem
            if rid not in t_rxn_ids:
                urxn_copy = urxn.copy()
                for met in urxn_copy.metabolites:
                    if met.id not in t_met_ids:
                        test_without.add_metabolites([met.copy()])
                        t_met_ids.add(met.id)
                test_without.add_reactions([urxn_copy])
                t_rxn_ids.add(rid)

        test_without.objective = 'biomass_human'
        tw_sol = test_without.optimize()
        tw_bm = tw_sol.objective_value if tw_sol.objective_value is not None else 0.0
        if tw_bm < 1e-6:
            needed_subsystems.add(ss)

    print(f"    Essential subsystems for biomass: {len(needed_subsystems)}")
    for ss in sorted(needed_subsystems):
        n_rxns = len(subsystem_rxns.get(ss, []))
        print(f"      {ss}: {n_rxns} reactions")

    # Now find minimal reaction set within needed subsystems
    # Add all needed subsystem reactions
    final_model = sm.copy()
    final_added = []
    f_met_ids = set(m.id for m in final_model.metabolites)
    f_rxn_ids = set(r.id for r in final_model.reactions)

    needed_rxns = []
    for ss in needed_subsystems:
        needed_rxns.extend(subsystem_rxns.get(ss, []))

    for rid in needed_rxns:
        if rid not in f_rxn_ids:
            urxn = universal.reactions.get_by_id(rid).copy()
            for met in urxn.metabolites:
                if met.id not in f_met_ids:
                    final_model.add_metabolites([met.copy()])
                    f_met_ids.add(met.id)
            final_model.add_reactions([urxn])
            f_rxn_ids.add(rid)
            final_added.append(rid)

    final_model.objective = 'biomass_human'
    final_sol = final_model.optimize()
    final_bm = final_sol.objective_value if final_sol.objective_value is not None else 0.0
    print(f"    Subsystem-based gap-fill: +{len(final_added)} rxns → biomass={final_bm:.6f}")

    if final_bm < 1e-6:
        # Subsystem approach didn't work perfectly, add ALL candidate reactions
        print(f"    Subsystem filtering insufficient, using all reference reactions...")
        final_model = test_model  # use the model with all candidates
        final_added = added_all if test_bm > 1e-6 else added_all2 if 'added_all2' in dir() else added_all
        final_bm = test_bm

    if final_bm < 1e-6:
        print(f"    *** Could not restore biomass for {sample_name}")
        gapfill_report.append({
            'sample': sample_name,
            'status': 'FAILED',
            'reactions_added': len(final_added),
            'biomass_after': final_bm,
            'details': 'Gap-fill unsuccessful'
        })
        continue

    # Try to reduce: remove reactions one by one if biomass is maintained
    print(f"    Pruning unnecessary reactions...")
    essential_gapfill = list(final_added)
    pruned = 0
    for rid in list(final_added):
        if rid not in set(r.id for r in final_model.reactions):
            essential_gapfill.remove(rid)
            continue
        # Try removing this reaction
        rxn_to_remove = final_model.reactions.get_by_id(rid)
        final_model.remove_reactions([rxn_to_remove])
        final_model.objective = 'biomass_human'
        prune_sol = final_model.optimize()
        prune_bm = prune_sol.objective_value if prune_sol.objective_value is not None else 0.0
        if prune_bm > 1e-6:
            # Reaction not needed
            essential_gapfill.remove(rid)
            pruned += 1
        else:
            # Reaction is essential, add it back
            urxn = universal.reactions.get_by_id(rid).copy()
            for met in urxn.metabolites:
                if met.id not in set(m.id for m in final_model.metabolites):
                    final_model.add_metabolites([met.copy()])
            final_model.add_reactions([urxn])

    final_model.objective = 'biomass_human'
    final_sol = final_model.optimize()
    final_bm = final_sol.objective_value if final_sol.objective_value is not None else 0.0

    print(f"    Pruned {pruned} unnecessary reactions")
    print(f"    Final: +{len(essential_gapfill)} essential gap-fill rxns → biomass={final_bm:.6f}")

    # Print the essential gap-fill reactions
    print(f"    Essential gap-fill reactions:")
    for rid in essential_gapfill:
        urxn = universal.reactions.get_by_id(rid)
        ss = urxn.subsystem if urxn.subsystem else "Unknown"
        print(f"      {rid}: {urxn.name} [{ss}]")

    # Save the gap-filled model
    out_path = os.path.join(MODELS_DIR, f"{sample_name}_model.xml")
    write_sbml_model(final_model, out_path)
    print(f"    Saved: {out_path}")

    gapfill_report.append({
        'sample': sample_name,
        'status': 'FIXED',
        'reactions_added': len(essential_gapfill),
        'biomass_after': final_bm,
        'details': "; ".join(essential_gapfill)
    })

# ── 6. Summary ──
print("\n" + "=" * 70)
print("[6] Gap-Fill Summary")
print("=" * 70)

# Save diagnosis
diag_df = pd.DataFrame(all_diagnoses)
diag_df.to_csv(os.path.join(RESULTS_DIR, "biomass_diagnosis.csv"), index=False)
print(f"  Diagnosis saved: biomass_diagnosis.csv")

# Save gap-fill report
report_df = pd.DataFrame(gapfill_report)
report_df.to_csv(os.path.join(RESULTS_DIR, "gapfill_report.csv"), index=False)
print(f"  Gap-fill report saved: gapfill_report.csv")

print(f"\n  Results:")
for row in gapfill_report:
    print(f"    {row['sample']}: {row['status']} (+{row['reactions_added']} rxns, biomass={row['biomass_after']:.6f})")

# ── 7. Final biomass verification for ALL models ──
print("\n" + "=" * 70)
print("[7] Final Biomass Verification (All Models)")
print("=" * 70)

all_samples = [
    "Lung_Healthy", "Lung_COVID19",
    "Heart_Healthy", "Heart_COVID19",
    "Liver_Healthy", "Liver_COVID19",
    "Kidney_Healthy", "Kidney_COVID19",
    "Intestine_Healthy", "Intestine_COVID19"
]

biomass_results = []
for sample_name in all_samples:
    xml_path = os.path.join(MODELS_DIR, f"{sample_name}_model.xml")
    if not os.path.exists(xml_path):
        print(f"  {sample_name}: model not found")
        continue
    sm = read_sbml_model(xml_path)
    if len(sm.reactions) == 0:
        print(f"  {sample_name}: empty model (0 reactions)")
        biomass_results.append({'sample': sample_name, 'n_reactions': 0, 'biomass_flux': 0.0})
        continue
    sm.objective = 'biomass_human'
    sol = sm.optimize()
    bm = sol.objective_value if sol.objective_value is not None else 0.0
    n_rxns = len(sm.reactions)
    status = "OK" if bm > 1e-6 else "ZERO"
    print(f"  {sample_name:25s}  rxns={n_rxns:5d}  biomass={bm:.6f}  [{status}]")
    biomass_results.append({'sample': sample_name, 'n_reactions': n_rxns, 'biomass_flux': bm})

# Save updated biomass fluxes
bm_df = pd.DataFrame(biomass_results)
bm_df.to_csv(os.path.join(RESULTS_DIR, "biomass_fluxes.csv"), index=False)
print(f"\n  Updated biomass_fluxes.csv saved")

total_time = time.time() - t0
print(f"\n  Total time: {total_time:.1f}s")
print("=" * 70)
