from setuptools import setup, find_packages

setup(
    name='troppo',
    version='0.2.0',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    install_requires=[
        "cobamp==0.2.1",
        "cobra==0.24.0",
        "numpy",
        "pandas",
        "matplotlib",
    ],
    extras_require={
        'full': [
            "scikit-learn",
            "scipy",
            "xlrd==1.2.0",
        ],
    },
    python_requires='>=3.8',
    author='Jorge Ferreira & Vítor Vieira',
    author_email='jorge.ferreira@ceb.uminho.pt',
    description='TROPPO - tINIT-based tissue-specific metabolic model reconstruction',
    license='GNU General Public License v3.0',
    keywords='metabolic model tINIT tissue-specific reconstruction omics',
    url='https://github.com/BioSystemsUM/troppo',
    long_description=open('README.rst').read(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
