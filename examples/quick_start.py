"""
Troppo Quick Start - Minimal tINIT Reconstruction Example

Prerequisites:
    - SBML metabolic model (e.g., Human-GEM.xml from https://github.com/SysBioChalmers/Human-GEM)
    - Expression data CSV with gene IDs as first column, samples as other columns
    - CPLEX or GUROBI solver installed and configured

Usage:
    python quick_start.py
"""

from troppo import TINITReconstructor

# ---- Step 1: Create reconstructor ----
recon = TINITReconstructor(
    model_path="path/to/Human-GEM.xml",       # Replace with your SBML model path
    expression_path="path/to/expression.csv",   # Replace with your expression CSV path
    solver="CPLEX",                             # or "GUROBI"
)

# ---- Step 2: Preview data before running ----
info = recon.preview_data()
print(f"Model: {info['model_reactions']} reactions, {info['model_genes']} genes")
print(f"Expression: {info.get('expression_genes', 'N/A')} genes, "
      f"{info.get('expression_samples', 'N/A')} samples")
print(f"Gene ID match rate: {info.get('gene_match_rate', 'N/A')}%")

# ---- Step 3: Run reconstruction ----
results = recon.run()

# ---- Step 4: View results ----
print("\n--- Summary ---")
print(results.summary())

print(f"\nCommon reactions across all samples: {len(results.get_common_reactions())}")

# ---- Step 5: Generate figures ----
from troppo.visualization import plot_overview
fig = plot_overview(results, save_dir="figures/")
print("\nOverview figure saved to figures/overview.pdf")
