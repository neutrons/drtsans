"""
    Setup script for ORNL SANS reduction
"""
from __future__ import absolute_import, division, print_function
from setuptools import setup, find_packages
import ornl


setup(name="ornl",
      version=ornl.__version__,
      description = "ORNL SANS reduction",
      url = "https://http://www.mantidproject.org",
      long_description = """ORNL SANS reduction""",
      license = "Apache License 2.0",
      zip_safe=False,
      packages=find_packages(),
      package_dir={},
      package_data={},
      install_requires=[],
      setup_requires=[],
)
