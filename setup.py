from setuptools import setup, find_packages
from assembly_finder import __version__

setup(
    name="assembly_finder",
    version=__version__,
    description="Snakemake pipeline for downloading assemblies from NCBI",
    url="https://github.com/metagenlab/assembly_finder",
    author="Farid Chaabane & Trestan Pillonel",
    author_email="trestan.pillonel@gmail.com",
    packages=find_packages(),
    package_data={"assembly_finder": ["Snakefile", "bin/*"]},
    data_files=[(".", ["LICENSE", "README.md"])],
    entry_points={
        "console_scripts": ["assembly_finder = assembly_finder.assembly_finder:cli"]
    },
)
