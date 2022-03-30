from setuptools import setup, find_packages
import glob
import os
import pkg_resources
from assembly_finder import __version__

setup(name='assembly_finder',
      version=__version__,
      packages=find_packages(),
      package_data={"":["assembly_finder/*", ]},
      description='Download assemblies from NCBI',
      url='https://github.com/metagenlab/assembly_finder',
      author='Farid Chaabane & Trestan Pillonel',
      author_email='trestan.pillonel@gmail.com',
      entry_points="""
      [console_scripts]
     assembly_finder = assembly_finder.assembly_finder:cli
      """,
      include_package_data=True,
      keywords=[],
      zip_safe=False)
