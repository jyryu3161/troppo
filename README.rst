
|License| |PyPI version|

TROPPO
============

*Troppo* is a Python package for **tINIT-based tissue-specific metabolic model reconstruction** using omics data.
It provides a streamlined pipeline from expression data to context-specific models with publication-quality analysis and visualization.

A (MI)LP solver (CPLEX or Gurobi) is required. The package uses optlang for solver abstraction.

Quick Start
~~~~~~~~~~~

::

    from troppo import TINITReconstructor

    # 3 lines: load, run, analyze
    recon = TINITReconstructor(
        model_path="Human-GEM.xml",
        expression_path="expression_data.csv",
        solver="CPLEX"
    )
    results = recon.run()
    results.summary()

Key Features
~~~~~~~~~~~~

**One-Stop tINIT Pipeline**

    - ``TINITReconstructor``: single entry point from raw data to reconstructed models
    - Automatic gene ID detection and conversion (Entrez, Ensembl, Symbol, HGNC)
    - Multi-sample batch processing (10-100 samples)
    - Configurable score integration strategies

**Rich Result Analysis**

    - ``ReconstructionResults``: container with built-in analysis methods
    - Reaction overlap, Jaccard similarity, subsystem coverage
    - Gene essentiality validation (accuracy, precision, recall, F1, MCC)
    - Theoretical yield calculations
    - Export to CSV, SBML, or pickle

**Publication-Quality Visualization**

    - 9 plot functions with Nature-style color palette, 300 DPI output
    - Model size comparison, reaction overlap heatmaps, subsystem coverage
    - PCA of model composition, biomass distribution, reaction frequency
    - ``plot_overview()``: composite 2x3 figure for quick reporting

**Flexible Gene ID Handling**

    - Automatic detection and conversion using HGNC database
    - Supports: ``entrez_id``, ``ensembl_gene_id``, ``symbol``, ``hgnc_id``, ``uniprot_ids``
    - Works with CSV/TSV expression data or pandas DataFrames

Installation
~~~~~~~~~~~~

From PyPI (stable)::

    pip install troppo

From GitHub (latest)::

    pip install git+https://github.com/BioSystemsUM/troppo.git

Dependencies::

    cobamp==0.2.1
    cobra==0.24.0
    numpy
    pandas
    matplotlib

Usage
~~~~~

**Basic Reconstruction**

::

    from troppo import TINITReconstructor

    recon = TINITReconstructor(
        model_path="Human-GEM.xml",
        expression_path="expression_data.csv",
        solver="CPLEX"
    )

    # Preview data before running
    recon.preview_data()

    # Configure tINIT parameters
    recon.configure(
        score_integration="adjusted",
        essential_reactions=["biomass_human"],
        and_func=min,
        or_func=max
    )

    # Run reconstruction
    results = recon.run()

**Analyzing Results**

::

    # Summary statistics
    results.summary()

    # Per-sample analysis
    active = results.get_active_reactions("sample1")
    common = results.get_common_reactions()
    unique = results.get_unique_reactions("sample1")

    # Cross-sample comparison
    results.jaccard_similarity()
    results.subsystem_coverage()
    results.biomass_fluxes()

    # Get a COBRA model for a specific sample
    model = results.get_model("sample1")

**Visualization**

::

    from troppo.visualization import (
        plot_model_sizes,
        plot_reaction_overlap_heatmap,
        plot_subsystem_coverage,
        plot_overview
    )

    # Individual plots
    fig = plot_reaction_overlap_heatmap(results)
    fig.savefig("jaccard_heatmap.pdf", dpi=300, bbox_inches="tight")

    # Composite overview (2x3 subplots)
    plot_overview(results, save_dir="figures/", format="pdf")

**Validation**

::

    # Gene essentiality
    ess_result = results.validate_gene_essentiality(
        essential_genes=["GAPDH", "PGK1"],
        non_essential_genes=["BRCA1", "TP53"]
    )

    # Theoretical yields
    yield_result = results.validate_theoretical_yields(
        substrate="glucose",
        product="biomass"
    )

**Working with Gene IDs**

Troppo automatically detects and converts gene IDs::

    from troppo.io import detect_gene_id_type, convert_gene_ids

    # Auto-detect ID type
    id_type = detect_gene_id_type(["2597", "3156", "5230"])
    # Returns: "entrez_id"

    # Convert between nomenclatures
    mapping = convert_gene_ids(
        ids=["2597", "3156"],
        from_type="entrez_id",
        to_type="symbol"
    )

When using ``TINITReconstructor``, set ``gene_id_type="auto"`` (default) to
let Troppo detect and convert IDs automatically::

    recon = TINITReconstructor(
        model_path="Human-GEM.xml",
        expression_path="expression_data.csv",
        gene_id_type="auto"  # auto-detect and convert
    )

**Saving and Loading Results**

::

    # Save results
    results.save("my_results.pkl")
    results.to_csv("reaction_matrix.csv")
    results.to_sbml("models/")

    # Load results
    from troppo import ReconstructionResults
    results = ReconstructionResults.load("my_results.pkl")

API Reference
~~~~~~~~~~~~~

**Core Classes**

    - ``troppo.TINITReconstructor`` - Main pipeline class
    - ``troppo.ReconstructionResults`` - Result container with analysis methods

**I/O Functions**

    - ``troppo.io.load_expression_csv()`` - Load CSV/TSV expression data
    - ``troppo.io.load_model()`` - Load SBML/JSON metabolic models
    - ``troppo.io.detect_gene_id_type()`` - Detect gene ID nomenclature
    - ``troppo.io.convert_gene_ids()`` - Convert between gene ID types

**Visualization Functions**

    - ``troppo.visualization.plot_model_sizes()`` - Active reactions per sample
    - ``troppo.visualization.plot_reaction_overlap_heatmap()`` - Jaccard similarity heatmap
    - ``troppo.visualization.plot_subsystem_coverage()`` - Subsystem activation heatmap
    - ``troppo.visualization.plot_biomass_distribution()`` - Biomass flux distribution
    - ``troppo.visualization.plot_pca_models()`` - PCA of reaction composition
    - ``troppo.visualization.plot_essentiality_metrics()`` - ROC curve and confusion matrix
    - ``troppo.visualization.plot_yield_comparison()`` - Predicted vs experimental yields
    - ``troppo.visualization.plot_reaction_frequency()`` - Reaction activation frequency
    - ``troppo.visualization.plot_overview()`` - Composite 2x3 report figure

Examples
~~~~~~~~

    - ``examples/quick_start.py`` - Minimal usage (5 lines)
    - ``examples/multi_sample_analysis.py`` - Multi-sample analysis with visualization

Credits and License
~~~~~~~~~~~~~~~~~~~

Developed at the Centre of Biological Engineering, University of Minho.

Released under the GNU Public License (version 3.0).


.. |License| image:: https://img.shields.io/badge/license-GPL%20v3.0-blue.svg
   :target: https://opensource.org/licenses/GPL-3.0
.. |PyPI version| image:: https://badge.fury.io/py/troppo.svg
   :target: https://badge.fury.io/py/troppo
