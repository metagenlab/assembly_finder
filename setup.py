from setuptools import setup, find_packages
from assembly_finder import __version__

setup(
    name="assembly_finder",
    version=__version__,
    description="Snakemake cli for downloading genomes with NCBI datasets",
    url="https://github.com/metagenlab/assembly_finder",
    python_requires=">=3.9",
    author="Farid Chaabane, Trestan Pillonel",
    author_email="farid.chaabane@chuv.ch, trestan.pillonel@gmail.com",
    packages=find_packages(),
    package_data={"assembly_finder": ["Snakefile", "bin/*"]},
    py_modules=["assembly_finder"],
    install_requires=[
        "rich-click>=1.7.4",
    ],
    data_files=[(".", ["LICENSE", "README.md"])],
    entry_points={
        "console_scripts": ["assembly_finder = assembly_finder.__main__:cli"]
    },
)
