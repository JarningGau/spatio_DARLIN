import os
import sys
from pathlib import Path

from setuptools import find_packages, setup

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))

setup(
    name="spatio-darlin",
    version="1.0.0",
    python_requires=">=3.8",
    packages=find_packages(),
    install_requires=[
        "darlin-core @ git+https://github.com/JarningGau/darlin-core.git",
        "numpy>=1.20.0",
        "pandas>=1.3.0",
    ],
    author="Shou-Wen Wang",
    author_email="wangshouwen@westlake.edu.cn",
    description="DARLIN snakemake pipeline for spatial lineage tracing",
    long_description=Path("README.md").read_text("utf-8"),
    license="BSD",
)
