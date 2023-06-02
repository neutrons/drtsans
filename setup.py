"""
Setup script for ORNL SANS reduction
"""
# from __future__ import absolute_import, division, print_function
import os
from setuptools import setup
from versioningit import get_cmdclasses

THIS_DIR = os.path.dirname(__file__)

# get list of scripts to install
scripts = [
    os.path.join(root, f)
    for root, _, files in os.walk(os.path.join(THIS_DIR, "scripts"))
    for f in files
    if f.endswith(".py")
]


def read_requirements_from_file(filepath):
    """Read a list of requirements from the given file and split into a
    list of strings. It is assumed that the file is a flat
    list with one requirement per line.
    :param filepath: Path to the file to read
    :return: A list of strings containing the requirements
    """
    with open(filepath, "r") as req_file:
        return req_file.readlines()


install_requires = read_requirements_from_file(os.path.join(THIS_DIR, "requirements.txt"))
test_requires = read_requirements_from_file(os.path.join(THIS_DIR, "requirements_dev.txt"))

setup(
    cmdclass=get_cmdclasses(),
    install_requires=install_requires,
)
