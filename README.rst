
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

**Benchmark Framework**

Systematically compare multiple reconstruction methods with comprehensive metrics and validation:

::

    from troppo.benchmark import BenchmarkRunner

    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'tinit', 'imat', 'fastcore'],
        biomass_reaction='biomass_human',
        essential_genes=['ENSG00000111640', 'ENSG00000100292'],  # GAPDH, HMGCR
        non_essential_genes=['ENSG00000139618', 'ENSG00000141510']  # BRCA2, TP53
    )

    comparison = runner.run_benchmark(validate_essentiality=True)

Features include:

    - Performance metrics (time, memory, model size)
    - Biological validation (biomass flux, gene essentiality)
    - Statistical metrics (accuracy, precision, recall, F1, MCC)
    - Rich visualizations (heatmaps, confusion matrices, radar charts)
    - Automated report generation

**Gene ID Support**

Troppo automatically handles multiple gene ID formats:

Supported ID types:

    - ``ensembl_gene_id`` - Ensembl Gene IDs (e.g., 'ENSG00000111640')
    - ``entrez_id`` - NCBI Entrez Gene IDs (e.g., '2597')
    - ``symbol`` - HGNC Gene Symbols (e.g., 'GAPDH')
    - ``hgnc_id`` - HGNC IDs (e.g., 'HGNC:4141')
    - ``uniprot_ids`` - UniProt IDs

Features:

    - Automatic ID type detection
    - Automatic conversion using HGNC database
    - Support for mixed ID types in workflows
    - Works seamlessly with benchmarking

**Quick Start: COVID-19 Case Study**

Here's a complete, runnable example using the COVID-19 case study data:

::

    import pandas as pd
    import cobra
    import re
    from troppo.methods_wrappers import ModelBasedWrapper
    from troppo.omics.readers.generic import TabularReader
    from troppo.benchmark import BenchmarkRunner

    # 1. Load model and expression data (Ensembl Gene IDs)
    model = cobra.io.read_sbml_model('examples/covid19_case_study/HumanGEM_Consistent_COVID19_HAM.xml')
    expression_data = pd.read_csv('examples/covid19_case_study/Desai-GTEx_ensembl.csv', index_col=0)

    # 2. Create data containers
    reader = TabularReader(
        path_or_df=expression_data,
        nomenclature='ensemble_gene_id',
        omics_type='transcriptomics'
    )
    omics_container = reader.to_containers()[0]

    # GPR parsing helper
    patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
    replace_alt_transcripts = lambda x: patt.sub('', x)

    # 3. Create model wrapper and map genes to reactions
    model_wrapper = ModelBasedWrapper(
        model=model,
        ttg_ratio=9999,
        gpr_gene_parse_function=replace_alt_transcripts
    )

    data_map = omics_container.get_integrated_data_map(
        model_reader=model_wrapper.model_reader,
        and_func=min,
        or_func=sum
    )

    # 4. Run benchmark comparing multiple methods
    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'fastcore', 'imat'],
        biomass_reaction='biomass_human'
    )

    comparison = runner.run_benchmark()

    # 5. View results
    print(comparison.get_summary_dataframe())

**Using tINIT with Metabolic Tasks**

For tINIT reconstruction, you can use metabolic task definitions:

::

    from troppo.methods.reconstruction.init import INIT, INITProperties
    import json

    # Load metabolic tasks (for tINIT)
    with open('examples/covid19_case_study/HumanGEM_nl2019_tasks_compact.json', 'r') as f:
        tasks = json.load(f)

    # Or use essential metabolic tasks from MCF7 study
    # tasks_file = 'examples/mcf7_case_study/data/metabolicTasks_Essential.xlsx'

    # Run tINIT with tasks
    properties = INITProperties(
        exp_vector=list(data_map.get_scores().values()),
        tasks=tasks
    )

    tinit = INIT(
        S=model_wrapper.S,
        lb=model_wrapper.lb,
        ub=model_wrapper.ub,
        properties=properties
    )

    result = tinit.run()

**Additional Examples and Documentation**

    - ``examples/covid19_case_study/`` - Complete COVID-19 case study
    - ``examples/mcf7_case_study/`` - MCF7 cell line case study
    - ``tests/Troppo_tutorial_omics_integration.ipynb`` - Comprehensive tutorial
    - ``tests/Troppo_tutorial_benchmark.ipynb`` - Benchmarking tutorial
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

**Basic Workflow**

1. **Load your data** - Use the COVID-19 example as a template::

    import pandas as pd
    import cobra
    from troppo.methods_wrappers import ModelBasedWrapper
    from troppo.omics.readers.generic import TabularReader

    # Load metabolic model and expression data
    model = cobra.io.read_sbml_model('path/to/model.xml')
    expression_data = pd.read_csv('path/to/expression.csv', index_col=0)

2. **Create data map** - Troppo automatically handles Ensembl, Entrez, Symbol, and other gene IDs::

    reader = TabularReader(
        path_or_df=expression_data,
        nomenclature='ensemble_gene_id',  # or 'entrez_id', 'symbol', etc.
        omics_type='transcriptomics'
    )
    omics_container = reader.to_containers()[0]
    model_wrapper = ModelBasedWrapper(model)
    data_map = omics_container.get_integrated_data_map(
        model_reader=model_wrapper.model_reader,
        and_func=min,
        or_func=sum
    )

3. **Run reconstruction** - Use any method (GIMME, tINIT, iMAT, FastCORE, CORDA, SWIFTCORE)::

    from troppo.benchmark import BenchmarkRunner

    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'fastcore', 'imat', 'tinit']
    )
    comparison = runner.run_benchmark()
    print(comparison.get_summary_dataframe())

**Quick Tips**

- All gene ID formats are supported (Ensembl, Entrez, Symbol, HGNC, UniProt)
- Use the COVID-19 case study in ``examples/covid19_case_study/`` as a starting point
- For tINIT, provide metabolic tasks (see example above)
- Compare methods using ``BenchmarkRunner`` to find the best one for your data

For detailed guides, see ``OMICS_INTEGRATION_GUIDE.md``

Credits and License
~~~~~~~~~~~~~~~~~~~

Developed at the Centre of Biological Engineering, University of Minho

Released under the GNU Public License (version 3.0).


.. |License| image:: https://img.shields.io/badge/license-GPL%20v3.0-blue.svg
   :target: https://opensource.org/licenses/GPL-3.0
.. |PyPI version| image:: https://badge.fury.io/py/troppo.svg
   :target: https://badge.fury.io/py/troppo
