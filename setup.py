"""
Setup script for ORNL SANS reduction
"""
from __future__ import absolute_import, division, print_function
import os
from setuptools import setup, find_packages
import versioneer

THIS_DIR = os.path.dirname(__file__)

# get list of scripts to install
scripts = [item for item in os.listdir(os.path.join(THIS_DIR, 'scripts'))
           if item.endswith('.py')]
# prepend directory
scripts = [os.path.join('scripts', item) for item in scripts]


def read_requirements_from_file(filepath):
    '''Read a list of requirements from the given file and split into a
    list of strings. It is assumed that the file is a flat
    list with one requirement per line.
    :param filepath: Path to the file to read
    :return: A list of strings containing the requirements
    '''
    with open(filepath, 'rU') as req_file:
        return req_file.readlines()


install_requires = read_requirements_from_file(os.path.join(THIS_DIR, 'requirements.txt'))
test_requires = read_requirements_from_file(os.path.join(THIS_DIR, 'requirements_dev.txt'))


setup(name="drtsans",
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description="Data Reduction Toolkit SANS reduction",
      url="https://http://www.mantidproject.org",
      long_description="""ORNL SANS reduction""",
      license="Apache License 2.0",
      zip_safe=False,
      packages=find_packages(),
      scripts=scripts,
      package_dir={},
      package_data={},
      install_requires=install_requires,
      setup_requires=['pytest-runner'],
      tests_require=test_requires,
      )
