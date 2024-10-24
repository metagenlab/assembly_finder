import os
from setuptools import setup, find_namespace_packages


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "assembly_finder",
            "assembly_finder.VERSION",
        )
    ) as f:
        return f.readline().strip()


def get_description():
    with open("README.md", "r") as fh:
        long_description = fh.read()
    return long_description


def get_data_files():
    data_files = [(".", ["README.md"])]
    return data_files


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="assembly_finder",
    packages=find_namespace_packages(),
    url="https://github.com/metagenlab/assembly_finder",
    python_requires=">=3.11",
    description="Snakemake-powered cli to download genomes using NCBI datasets",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Farid Chaabane & Trestan Pillonel",
    author_email="farid.chaabane@chuv.ch, trestan.pillonel@gmail.com",
    data_files=get_data_files(),
    py_modules=["assembly_finder"],
    install_requires=[
        "snakemake>=8.0.0",
        "snaketool-utils>=0.0.5",
        "attrmap>=0.0.7",
        "pyyaml>=6.0",
        "pandas>=2.2.1",
        "rich-click>=1.8.3",
    ],
    entry_points={"console_scripts": ["assembly_finder=assembly_finder.__main__:main"]},
    include_package_data=True,
)
