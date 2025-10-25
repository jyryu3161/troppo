
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

    from troppo.benchmark import quick_benchmark

    # Compare methods with one line
    summary = quick_benchmark(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'tinit', 'imat', 'fastcore']
    )

Features include:

    - Automatic performance metrics (time, memory, model size)
    - Biological validation (biomass flux, task completion)
    - Network quality metrics (consistency, blocked reactions)
    - Rich visualizations (heatmaps, radar charts, Pareto fronts)
    - Automated report generation

**Tutorials and Examples**

    - ``tests/Troppo_tutorial_omics_integration.ipynb`` - Comprehensive omics integration tutorial
    - ``tests/Troppo_tutorial_benchmark.ipynb`` - Method comparison and benchmarking
    - ``examples/custom_method_example.py`` - Complete custom method implementation
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
