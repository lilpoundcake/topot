#!/usr/bin/env python
"""Setup script for topot package"""

from setuptools import setup, find_packages
from pathlib import Path

# Read long description from README
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="topot",
    version="0.1.0",
    description="Extract lambda-specific structures from GROMACS dual topology files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Contributors",
    author_email="topot@example.com",
    url="https://github.com/yourusername/topot",
    license="MIT",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8.1",
    install_requires=[
        "numpy>=1.19",
        "biopython>=1.81",
    ],
    extras_require={
        "dev": [
            "pytest>=7.4",
            "pytest-cov>=4.0",
            "black>=23.0",
            "flake8>=6.0",
            "mypy>=1.5",
            "isort>=5.0",
        ],
        "docs": [
            "sphinx>=5.0",
            "sphinx-rtd-theme>=1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "topot=topot.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="GROMACS molecular-dynamics free-energy PMX topology protein-mutation",
)
