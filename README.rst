
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

Key Features
~~~~~~~~~~~~

**Flexible Gene ID Handling**

    - Automatic detection and conversion of gene IDs (Entrez, Ensembl, Symbol, HGNC, UniProt)
    - No manual ID conversion required
    - Works with any expression data format
    - Seamless integration with metabolic models

**Comprehensive Benchmarking**

    - Compare multiple reconstruction methods
    - Performance metrics (time, memory, model size)
    - Biological validation (gene essentiality, theoretical yields)
    - Statistical metrics (accuracy, precision, recall, F1, MCC)
    - Rich visualizations and automated reports

**Extensible Architecture**

    - Plugin-like system for adding new methods
    - Easy integration with existing tools (COBRApy, FRAMED, CAMEO)
    - Modular design for custom workflows

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

**Expression Data with Any Gene ID Nomenclature**

Troppo now supports expression data with any gene ID format, automatically converting to match your model:

::

    from troppo.omics import create_data_map_from_dict

    # Expression data with Entrez IDs (or any other ID type)
    expression_data = {
        '2597': 100.5,   # GAPDH (Entrez ID)
        '3156': 85.2,    # HMGCR (Entrez ID)
        '5230': 95.8     # PGK1 (Entrez ID)
    }

    # Automatically converts IDs to match your model
    data_map = create_data_map_from_dict(
        expression_data,
        model_wrapper,
        gene_id_type='entrez_id',  # or auto-detect if omitted
        auto_convert=True,
        verbose=True
    )

**Multi-Nomenclature Gene ID Support**

Use different ID types for expression data and validation data in the same workflow:

::

    from troppo.omics import create_data_map_from_dict
    from troppo.benchmark import BenchmarkRunner

    # Expression data: Entrez IDs
    expression_data = {'2597': 100.5, '3156': 85.2}
    data_map = create_data_map_from_dict(
        expression_data,
        model_wrapper,
        gene_id_type='entrez_id',
        auto_convert=True
    )

    # Benchmark with different ID types for validation
    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'tinit', 'imat'],
        # Essential genes in Ensembl ID format
        essential_genes=['ENSG00000111640', 'ENSG00000102144'],
        essential_genes_id_type='ensembl_gene_id',
        # Non-essential genes in Symbol format
        non_essential_genes=['HMGCR', 'TP53'],
        non_essential_genes_id_type='symbol'
    )

    # All IDs automatically converted to match model
    comparison = runner.run_benchmark(validate_essentiality=True)

Supported ID types:

    - ``entrez_id`` - NCBI Entrez Gene IDs (e.g., '2597')
    - ``ensembl_gene_id`` - Ensembl Gene IDs (e.g., 'ENSG00000111640')
    - ``symbol`` - HGNC Gene Symbols (e.g., 'GAPDH')
    - ``hgnc_id`` - HGNC IDs (e.g., 'HGNC:4141')
    - ``uniprot_ids`` - UniProt IDs

Features include:

    - Automatic ID type detection for expression data
    - Automatic conversion using HGNC database
    - Support for mixed ID types in single workflow
    - Works seamlessly with benchmark and validation
    - Verbose logging of conversion process

**Quick Start**

Basic workflow for omics integration:

::

    from troppo.methods_wrappers import ModelBasedWrapper
    from troppo.omics import create_data_map_from_dict
    from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
    import cobra

    # 1. Load metabolic model
    model = cobra.io.read_sbml_model('path/to/model.xml')
    model_wrapper = ModelBasedWrapper(model)

    # 2. Prepare expression data (any gene ID format)
    expression_data = {
        '2597': 100.5,   # GAPDH
        '3156': 85.2,    # HMGCR
        '5230': 95.8     # PGK1
        # ... more genes
    }

    # 3. Create data map (automatic ID conversion)
    data_map = create_data_map_from_dict(
        expression_data,
        model_wrapper,
        gene_id_type='entrez_id',
        auto_convert=True
    )

    # 4. Run tissue-specific reconstruction
    properties = GIMMEProperties(
        exp_vector=list(data_map.get_scores().values()),
        obj_frac=0.9
    )

    gimme = GIMME(
        S=model_wrapper.S,
        lb=model_wrapper.lb,
        ub=model_wrapper.ub,
        properties=properties
    )

    result = gimme.run()

**Tutorials and Examples**

    - ``tests/Troppo_tutorial_omics_integration.ipynb`` - Comprehensive omics integration tutorial
    - ``tests/Troppo_tutorial_benchmark.ipynb`` - Method comparison and benchmarking
    - ``examples/custom_method_example.py`` - Complete custom method implementation
    - ``examples/benchmark_with_validation_example.py`` - Gene essentiality and yield validation
    - ``examples/benchmark_with_entrez_ids_example.py`` - Using Entrez IDs in benchmarks
    - ``examples/expression_data_with_entrez_ids.py`` - Expression data with different ID types
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

Usage
~~~~~

**Using Expression Data with Different Gene IDs**

Troppo makes it easy to work with expression data in any gene ID format:

Method 1 - Direct from dictionary (simplest)::

    from troppo.omics import create_data_map_from_dict

    expression_data = {'2597': 100.5, '3156': 85.2}  # Entrez IDs
    data_map = create_data_map_from_dict(
        expression_data,
        model_wrapper,
        gene_id_type='entrez_id',
        auto_convert=True
    )

Method 2 - Using OmicsContainer::

    from troppo.omics import OmicsContainer, create_compatible_data_map

    omics_container = OmicsContainer(
        omicstype='transcriptomics',
        condition='sample1',
        data=expression_data,
        nomenclature='entrez_id'
    )

    data_map = create_compatible_data_map(
        omics_container,
        model_wrapper,
        auto_convert=True
    )

**Running Benchmarks**

Compare multiple methods with automatic validation::

    from troppo.benchmark import BenchmarkRunner

    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'fastcore', 'imat'],
        essential_genes=['GAPDH', 'PGK1'],
        carbon_sources=['glucose', 'acetate']
    )

    comparison = runner.run_benchmark(
        validate_essentiality=True,
        validate_yields=True
    )

    # View results
    summary = comparison.get_summary_dataframe()
    print(summary)

**Troubleshooting Common Issues**

ID mismatch between expression data and model?
    Use ``auto_convert=True`` when creating data map

Want to check what ID type your data uses?::

    from troppo.omics import detect_expression_data_id_type
    id_type = detect_expression_data_id_type(expression_data)

Need to convert IDs manually?::

    omics_container.convertIds('symbol')  # Convert to gene symbols

For more details, see ``OMICS_INTEGRATION_GUIDE.md``

Credits and License
~~~~~~~~~~~~~~~~~~~

Developed at the Centre of Biological Engineering, University of Minho

Released under the GNU Public License (version 3.0).


.. |License| image:: https://img.shields.io/badge/license-GPL%20v3.0-blue.svg
   :target: https://opensource.org/licenses/GPL-3.0
.. |PyPI version| image:: https://badge.fury.io/py/troppo.svg
   :target: https://badge.fury.io/py/troppo
