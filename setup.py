#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2019--, CMI development team.

# ----------------------------------------------------------------------------
import re
import ast
from glob import glob
from setuptools import find_packages, setup

classes = """
    Development Status :: 3 - Alpha
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split("\n") if s]

description = "Downloading and formatting of data for Qiita"

with open("README.md") as f:
    long_description = f.read()

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r"__version__\s+=\s+(.*)")
with open("qebil/__init__.py", "rb") as f:
    hit = _version_re.search(f.read().decode("utf-8")).group(1)
    version = str(ast.literal_eval(hit))

standalone = ["qebil=qebil.commands:cli"]

setup(
    name="qebil",
    version=version,
    license="TBD",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Center for Microbiome Innovation",
    author_email="cmiinfo@ucsd.edu",
    maintainer="Center for Microbiome Innovation",
    maintainer_email="cmiinfo@ucsd.edu",
    packages=find_packages(),
    install_requires=[
        "pandas >= 0.23.4",
        "numpy >= 1.15.4",
        'requests',
        'PyPDF2',
        'click',
        'black',
        'xmltodict',
        'bs4',
        'pyyaml >=5.1',
        'xlrd == 1.2.0'
    ],
    classifiers=classifiers,
    entry_points={"console_scripts": standalone},
    include_package_data=True,
    package_data={'qebil': ['support_files/validators/*','support_files/*']},
    zip_safe=False,
)
