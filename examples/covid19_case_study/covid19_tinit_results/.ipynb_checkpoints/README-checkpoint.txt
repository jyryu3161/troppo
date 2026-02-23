COVID-19 Metabolic Model Reconstruction Results
===============================================
Date: 2026-02-20
Status: Pending re-run with biomass chain fix

1. Project Overview
-------------------
- Method: tINIT (Task-driven Integrative Network Inference for Tissues)
- Input: HumanGEM (Universal), COVID-19 vs Healthy RNA-seq (Desai et al.)
- Output: 10 Context-Specific Metabolic Models (5 Healthy, 5 COVID-19)
- Tissues: Lung, Heart, Liver, Kidney, Intestine
- Thresholding: Local T2 strategy (quantiles: 0.50 lower, 0.75 upper, 0.50 local)

2. Reconstruction Metrics (Previous Run)
-----------------------------------------
- Total reactions in universal model: 10,348
- Reaction retention per model:
    Lung_Healthy       6,418 rxns  (62.0%)
    Lung_COVID19       6,378 rxns  (61.6%)
    Heart_Healthy      6,021 rxns  (58.2%)
    Heart_COVID19      5,848 rxns  (56.5%)
    Liver_Healthy      6,509 rxns  (62.9%)
    Liver_COVID19      6,730 rxns  (65.0%)
    Kidney_Healthy     6,315 rxns  (61.0%)
    Kidney_COVID19     5,835 rxns  (56.4%)
    Intestine_Healthy  6,571 rxns  (63.5%)
    Intestine_COVID19  6,301 rxns  (60.9%)

- Jaccard Similarity (Healthy vs COVID-19 within tissue):
    Lung:       0.869  (highest - robust core network)
    Intestine:  0.848
    Heart:      0.828
    Liver:      0.749
    Kidney:     0.634  (lowest - severe metabolic perturbation)

3. Biomass Flux Issue & Fix Applied
------------------------------------
Problem:
  All 10 models had biomass flux = 0.0 because tINIT removed critical
  biomass chain reactions during reconstruction.

Root Cause (3 layers):
  Layer 1: HMR_10023 + HMR_10024 (biomass transport/exchange) removed,
           causing temp001c (biomass product) to be a dead-end metabolite.
  Layer 2: HMR_8617 (triphosphate hydrolysis) removed, blocking
           cobamide-coenzyme via m03057m dead-end.
  Layer 3: HMR_10065 (cofactor pool assembly) required for biomass but
           removed in some models.
  Layer 4: BH4/THF/CL pathway gaps in 6/10 models (biopterin, folate
           metabolism missing reactions).

Fix Applied (test_covid19_tinit.py):
  1. Added HMR_10023, HMR_10024, HMR_8617, HMR_10065 to essential_reactions
     so tINIT never removes them.
  2. Added post-hoc gap-fill step using 16 reactions from
     diagnostic_universal_fix.csv (BH4/THF/CL pathways) for any model
     that still has zero biomass after reconstruction.

Status: Awaiting re-run to verify fix produces non-zero biomass.
  Command: python test_covid19_tinit.py

4. Biological Insights (from Network Topology)
-----------------------------------------------
Key metabolic shifts observed in COVID-19 vs Healthy:

A. Loss of Function (Downregulated in COVID-19):
   - Metabolism of other amino acids (Delta: -0.40)
   - Vitamin C metabolism (Delta: -0.40): impaired antioxidant defense
   - Biotin metabolism (Delta: -0.35): cofactor depletion
   - Vitamin B6 metabolism (Delta: -0.31): cofactor depletion
   - Cholesterol biosynthesis 2 (Delta: -0.30): altered sterol metabolism

B. Gain of Function (Upregulated in COVID-19):
   - Beta oxidation of odd-chain fatty acids, peroxisomal (Delta: +0.20)

C. Tissue-Specific Observations:
   - Kidney shows largest divergence (Jaccard 0.634), suggesting
     COVID-19 causes severe renal metabolic reprogramming.
   - Lung shows highest conservation (Jaccard 0.869), indicating
     a robust core metabolic network despite infection.
   - Liver_COVID19 has the most reactions (6,730), possibly reflecting
     hepatic metabolic compensation during infection.

5. Validation (Previous Run, Affected by Biomass=0)
-----------------------------------------------------
- Theoretical Yield: All 0.0 (expected, due to biomass chain issue)
- Gene Essentiality: Skipped for all models (WT growth too low)
- These will produce meaningful results after biomass fix re-run.

6. File Manifest
-----------------
Data Files:
  reconstruction_summary.csv    Model sizes and retention percentages
  reaction_matrix.csv           Binary reaction presence/absence (10 x 10,348)
  jaccard_similarity.csv        Pairwise Jaccard similarity (10 x 10)
  subsystem_coverage.csv        Metabolic pathway coverage per model
  differential_reactions.csv    Healthy-only / COVID-only / shared reactions
  tissue_comparison.csv         Per-tissue Healthy vs COVID-19 statistics
  biomass_fluxes.csv            Biomass flux per model (pending re-run)
  benchmark_yields.csv          Theoretical yield results
  benchmark_essentiality.csv    Gene essentiality validation
  reconstruction_results.pkl    Pickled ReconstructionResults object

Models:
  models/                       SBML/XML files for all 10 reconstructed models

Figures:
  figures/                      Visualization plots (PCA, heatmaps, etc.)

Supporting:
  ../diagnostic_universal_fix.csv   Gap-fill reaction list (16 rxns)
