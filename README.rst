.. image:: https://img.shields.io/readthedocs/drtsans.svg?logo=read-the-docs
   :target: https://drtsans.readthedocs.io/en/latest
   :alt: Documentation Status

.. image:: https://github.com/neutrons/drtsans/actions/workflows/package.yml/badge.svg?branch=main
   :alt: CI
   :target: https://github.com/neutrons/drtsans/actions/workflows/package.yml?query=branch:main

.. image:: https://codecov.io/gh/neutrons/drtsans/branch/next/graph/badge.svg?token=q1f07RUI88
   :alt: codecov
   :target: https://codecov.io/gh/neutrons/drtsans

.. image:: https://www.bestpractices.dev/projects/6400/badge
   :alt: CII Best Practices
   :target: https://www.bestpractices.dev/projects/6400



===========
drt-sans
===========

Data Reduction Toolkit for Small Angle Neutron Scattering
--------------------------------------------------------

This package is a collection of functionality for reducing SANS data developed in collaboration with the instrument
scientists at the High Flux Isotope Reactor (HFIR) and Spallation Neutron Source (SNS) at Oak Ridge National Laboratory.

While much of the functionality is generic, this implementation is aimed at reducing data from BIOSANS, EQSANS,
and GPSANS. As appropriate, this work is an abstraction layer on top of the Mantid project.

**This is a python3-only package.**

--------------------------
General User Documentation
--------------------------

Most of the current documentation here refers to the `drtsans` Python programming application interface (API), suitable for
power users. For the general user of `drtsans`, https://sans.ornl.gov is available. For developers and users wanting to have a more in depth knowledge about certain subjects, visit https://drtsans.readthedocs.io/latest/index.html.

The general user documentation is hosted in repository `https://code-int.ornl.gov/sites/sans`.

Visit the repository at `https://github.com/neutrons/drtsans`.

------------------------------
Usage from Provided Front-ends
------------------------------

For end users, go to the `QA version <http://scse-ui.ornl.gov:8080/>`_.

Use `Jupyter <https://jupyter.sns.gov/>`_ to have a play with the code.
The kernel to select is ``sans at ...``.

One can run scripts directly on the `analysis <https://analysis.sns.gov/>`_ cluster.
To do that, open a terminal and activate the desired conda environment. The options are:

* ``sans``: the latest stable release
* ``sans-qa``: the future stable release (to be tested right before the next iteration)
* ``sans-dev``: the latest development version

The easiest way to start an interactive IPython session is by running:

.. code-block:: shell

   $ drtsans

Adding ``--qa`` or ``--dev`` will start the QA or development version respectively.
The ``drtsans`` wrapper script launches IPython with the selected conda environment located in ``/opt/anaconda/envs/``
and deactivates the conda environment when the session ends.

One must have an XCAMS account to use either the Jupyter kernel provided above.

-----------------------------------------------
Set-up for Development in a Virtual Environment
-----------------------------------------------
