from setuptools import setup, find_packages
import os


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "assembly_finder",
            "assembly_finder.VERSION",
        )
    ) as f:
        return f.readline().strip()


setup(
    name="assembly_finder",
    version=get_version(),
    description="Snakemake cli for downloading genomes with NCBI datasets",
    url="https://github.com/metagenlab/assembly_finder",
    python_requires=">=3.9",
    author="Farid Chaabane, Trestan Pillonel",
    author_email="farid.chaabane@chuv.ch, trestan.pillonel@gmail.com",
    packages=find_packages(),
    package_data={"assembly_finder": ["Snakefile", "bin/*", "assembly_finder.VERSION"]},
    py_modules=["assembly_finder"],
    install_requires=[
        "rich-click>=1.7.4",
    ],
    data_files=[(".", ["LICENSE", "README.md"])],
    entry_points={
        "console_scripts": ["assembly_finder = assembly_finder.__main__:cli"]
    },
)
