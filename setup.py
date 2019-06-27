"""
    Setup script for ORNL SANS reduction
"""
from __future__ import absolute_import, division, print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import versioneer


class PyTest(TestCommand):
    user_options = [("pytest-args=", "a", "Arguments to pass to pytest")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ""

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        import shlex
        import sys

        errno = pytest.main(shlex.split(self.pytest_args))
        sys.exit(errno)


commands = versioneer.get_cmdclass()
commands['pytest'] = PyTest  # add the custom pytest command

setup(name="ornl",
      version=versioneer.get_version(),
      cmdclass=commands,
      description="ORNL SANS reduction",
      url="https://http://www.mantidproject.org",
      long_description="""ORNL SANS reduction""",
      license="Apache License 2.0",
      zip_safe=False,
      packages=find_packages(),
      package_dir={},
      package_data={},
      install_requires=[],
      setup_requires=[],
      )
