
|License| |PyPI version|

TROPPO
============

*Troppo* (Tissue-specific RecOnstruction and Phenotype Prediction using Omics data) is a Python package containing methods
for tissue specific reconstruction to use with constraint-based models. The main purpose of this package is to provide
an open-source framework which is both modular and flexible to be integrated with other packages, such as cobrapy, framed
or cameo whose already provide generic data structures to interpret metabolic models.

A (MI)LP solver is required to use most of the present methods. The current methods support optlang, which in turn allow
the use of solvers like CPLEX or Gurobi.

*Troppo*'s documentation is available at http://troppo-bisbi.readthedocs.io/.

The current methods implemented are:
    - FastCORE
    - CORDA
    - GIMME
    - (t)INIT
    - iMAT
    - SWIFTCORE

Methods to be implemented later:
    - MBA
    - mCADRE
    - PRIME

New Features
~~~~~~~~~~~~

**Extensible Plugin System**

Troppo now includes a plugin-like registry system that makes it easy to add new reconstruction methods:

::

    from troppo.methods.registry import MethodRegistry, register_method

    # Register your custom method
    @register_method('my_method', description='My custom method')
    class MyMethod(ContextSpecificModelReconstructionAlgorithm):
        # ... implementation ...

    # Use it immediately
    method = MethodRegistry.create_method('my_method', S, lb, ub, properties)
    result = method.run()

**Benchmark Framework**

Compare multiple methods systematically with comprehensive metrics:

::

    from troppo.benchmark import BenchmarkRunner

    # Compare methods with validation
    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'tinit', 'imat', 'fastcore'],
        # Gene essentiality validation
        essential_genes=['GAPDH', 'HMGCR', 'PGK1'],
        non_essential_genes=['BRCA1', 'BRCA2', 'TP53'],
        # Theoretical yield validation
        carbon_sources=['glucose', 'acetate', 'glycerol'],
        biomass_reaction='biomass_human'
    )

    comparison = runner.run_benchmark(
        validate_essentiality=True,
        validate_yields=True
    )

Features include:

    - Automatic performance metrics (time, memory, model size)
    - Biological validation (biomass flux, task completion)
    - Gene essentiality validation (accuracy, precision, recall, F1, MCC)
    - Theoretical yield calculations (aerobic and anaerobic conditions)
    - Network quality metrics (consistency, blocked reactions)
    - Rich visualizations (heatmaps, radar charts, Pareto fronts, confusion matrices)
    - Automated report generation

**Multi-Nomenclature Gene ID Support**

Troppo automatically handles different gene ID nomenclatures for validation:

::

    from troppo.benchmark import BenchmarkRunner

    # Use Entrez IDs for validation data
    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'tinit', 'imat'],
        # Essential genes in Entrez ID format
        essential_genes=['2597', '3156', '5230'],  # GAPDH, HMGCR, PGK1
        essential_genes_id_type='entrez_id',
        # Non-essential genes in Ensembl ID format
        non_essential_genes=['ENSG00000139618', 'ENSG00000141510'],
        non_essential_genes_id_type='ensembl_gene_id'
    )

    # IDs are automatically converted to match your model's nomenclature
    comparison = runner.run_benchmark(validate_essentiality=True)

Supported ID types:

    - ``entrez_id`` - NCBI Entrez Gene IDs (e.g., '2597')
    - ``ensembl_gene_id`` - Ensembl Gene IDs (e.g., 'ENSG00000111640')
    - ``symbol`` - HGNC Gene Symbols (e.g., 'GAPDH')
    - ``hgnc_id`` - HGNC IDs (e.g., 'HGNC:4141')
    - ``uniprot_ids`` - UniProt IDs

Features include:

    - Automatic ID type detection
    - Automatic conversion using HGNC database
    - Support for mixed ID types
    - Flexible ID handler for complex scenarios
    - Verbose logging of conversion process

**Tutorials and Examples**

    - ``tests/Troppo_tutorial_omics_integration.ipynb`` - Comprehensive omics integration tutorial
    - ``tests/Troppo_tutorial_benchmark.ipynb`` - Method comparison and benchmarking
    - ``examples/custom_method_example.py`` - Complete custom method implementation
    - ``examples/benchmark_with_validation_example.py`` - Gene essentiality and yield validation
    - ``examples/benchmark_with_entrez_ids_example.py`` - Using Entrez IDs in benchmarks
    - ``run_omics_integration.sh`` - Automated pipeline script

**Documentation**

    - ``OMICS_INTEGRATION_GUIDE.md`` - User guide for omics integration
    - ``EXTENSIBILITY_GUIDE.md`` - Developer guide for extending Troppo

Instalation from PyPI (stable releases)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    pip install troppo

Instalation from github (latest development release)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    pip install git+https://github.com/BioSystemsUM/troppo.git



Credits and License
~~~~~~~~~~~~~~~~~~~

Developed at the Centre of Biological Engineering, University of Minho

Released under the GNU Public License (version 3.0).


.. |License| image:: https://img.shields.io/badge/license-GPL%20v3.0-blue.svg
   :target: https://opensource.org/licenses/GPL-3.0
.. |PyPI version| image:: https://badge.fury.io/py/troppo.svg
   :target: https://badge.fury.io/py/troppo
