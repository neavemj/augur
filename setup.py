<<<<<<< .merge_file_unrNfV
import os
from setuptools import setup

setup(
        name = "augur",
        version = "0.1.0",
        author = "nextstrain developers",
        author_email = "trevor@bedford.io, richard.neher@unibas.ch",
        description = ("Pipelines for real-time phylogenetic analysis"),
        license = "MIT",
        keywords = "nextstrain, molecular epidemiology",
        url = "https://github.com/nextstrain/augur",
        packages=['augur'],
        install_requires = [
            "biopython >=1.69, ==1.*",
            "boto >=2.38, ==2.*",
            "cvxopt >=1.1.8, ==1.1.*",
            "ipdb >=0.10.1, ==0.10.*",
            "matplotlib >=2.0, ==2.*",
            "pandas >=0.16.2, <0.18.0",
            "pytest >=3.2.1, ==3.*",
            "seaborn >=0.6.0, ==0.6.*",
            "tox >=2.8.2, ==2.*",
            "treetime ==0.4.0"
        ],
        dependency_links = [
            "https://api.github.com/repos/neherlab/treetime/tarball/v0.4.0#egg=treetime-0.4.0"
        ],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Science",
            "License :: OSI Approved :: MIT License",
            ],
        scripts=['bin/augur']
        )
=======
"""
Pipeline components for real-time virus analysis
Author: Richard Neher and Trevor Bedford
"""
import os
from setuptools import setup, find_packages

setup(
    name = "augur",
    version = "1.0.0",
    author = "Richard Neher and Trevor Bedford",
    author_email = "trevor@bedford.io",
    description = ("Pipeline components for real-time virus analysis"),
    license = "GNU Affero General Public License v3.0",
    keywords = "",
    url = "https://github.com/nextstrain/augur",
    packages=find_packages(where='src'),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Science",
        "License :: GNU Affero General Public License v3.0",
    ]
)
>>>>>>> .merge_file_fUdj8V
