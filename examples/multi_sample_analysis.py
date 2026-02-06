"""
Multi-Sample Analysis with tINIT

Demonstrates:
    1. Loading expression data with automatic gene ID conversion
    2. Configuring reconstruction parameters
    3. Running reconstruction for multiple samples
    4. Analyzing results (Jaccard similarity, subsystem coverage, PCA)
    5. Generating individual publication-quality figures
    6. Exporting results

Prerequisites:
    - SBML metabolic model
    - Expression data CSV (multiple samples)
    - CPLEX or GUROBI solver
    - Optional: scikit-learn for PCA plot
"""

import os
from troppo import TINITReconstructor
from troppo.visualization import (
    plot_model_sizes,
    plot_reaction_overlap_heatmap,
    plot_subsystem_coverage,
    plot_pca_models,
    plot_reaction_frequency,
    plot_overview,
)


def main():
    # ---- Configuration ----
    MODEL_PATH = "path/to/Human-GEM.xml"
    EXPRESSION_PATH = "path/to/expression_data.csv"
    SOLVER = "CPLEX"
    OUTPUT_DIR = "results"
    FIGURE_DIR = os.path.join(OUTPUT_DIR, "figures")

    os.makedirs(FIGURE_DIR, exist_ok=True)

    # ---- 1. Initialize reconstructor ----
    print("=" * 60)
    print("TROPPO - Multi-Sample tINIT Reconstruction")
    print("=" * 60)

    recon = TINITReconstructor(
        model_path=MODEL_PATH,
        expression_path=EXPRESSION_PATH,
        solver=SOLVER,
        verbose=True,
    )

    # ---- 2. Configure parameters ----
    recon.configure(
        production_weight=0.5,
        no_reverse_loops=True,
        score_integration="adjusted",
    )

    # Preview data
    info = recon.preview_data()
    print(f"\nData Preview:")
    for k, v in info.items():
        if k != 'sample_names':
            print(f"  {k}: {v}")

    # ---- 3. Run reconstruction ----
    results = recon.run()

    # ---- 4. Analyze results ----
    print("\n" + "=" * 60)
    print("ANALYSIS")
    print("=" * 60)

    # Summary
    summary = results.summary()
    print("\nPer-sample summary:")
    print(summary.to_string())

    # Common vs unique reactions
    common = results.get_common_reactions()
    print(f"\nReactions common to all samples: {len(common)}")

    for sample in results.sample_names:
        unique = results.get_unique_reactions(sample)
        print(f"  Reactions unique to {sample}: {len(unique)}")

    # Jaccard similarity
    jac = results.jaccard_similarity()
    print("\nJaccard similarity matrix:")
    print(jac.round(3).to_string())

    # ---- 5. Generate figures ----
    print(f"\nSaving figures to {FIGURE_DIR}/...")

    fig = plot_model_sizes(results, save_path=os.path.join(FIGURE_DIR, "model_sizes.pdf"))
    print("  model_sizes.pdf")

    fig = plot_reaction_overlap_heatmap(results, save_path=os.path.join(FIGURE_DIR, "jaccard_heatmap.pdf"))
    print("  jaccard_heatmap.pdf")

    fig = plot_subsystem_coverage(results, top_n=20, save_path=os.path.join(FIGURE_DIR, "subsystem_coverage.pdf"))
    print("  subsystem_coverage.pdf")

    fig = plot_reaction_frequency(results, top_n=30, save_path=os.path.join(FIGURE_DIR, "reaction_frequency.pdf"))
    print("  reaction_frequency.pdf")

    try:
        fig = plot_pca_models(results, save_path=os.path.join(FIGURE_DIR, "pca_models.pdf"))
        print("  pca_models.pdf")
    except ImportError:
        print("  pca_models.pdf (skipped - install scikit-learn)")

    fig = plot_overview(results, save_dir=FIGURE_DIR)
    print("  overview.pdf")

    # ---- 6. Export results ----
    results.to_csv(os.path.join(OUTPUT_DIR, "reaction_matrix.csv"))
    print(f"\nReaction matrix saved to {OUTPUT_DIR}/reaction_matrix.csv")

    results.save(os.path.join(OUTPUT_DIR, "results.pkl"))
    print(f"Results object saved to {OUTPUT_DIR}/results.pkl")

    print("\nDone!")


if __name__ == "__main__":
    main()
