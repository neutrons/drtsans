.. image:: https://www.bestpractices.dev/projects/6400/badge
   :target: https://www.bestpractices.dev/projects/6400

drt-sans
========

Data Reduction Toolkit for Small Angle Neutron Scattering

This packages is a collection of functionality for reducing SANS data developed in collaboration with the instrument
scientists at the High Flux Isotope Reactor (HFIR) and Spallation Neutron Source (SNS) at Oak Ridge National Laboratory.

While much of the functionality is generic, this implementation is aimed at reducing data from BIOSANS, EQSANS,
and GPSANS. As appropriate, this work is an abstraction layer on top of the mantid project.

**This is a python3 only package.**

General User documentation
--------------------------

Most of the current documentation here refers to `drtsans` python programming application interface (API), suitable for
power users. For the general user of `drtsans`, https://sans.ornl.gov is available.

The general user documentation is hosted in repository https://code-int.ornl.gov/sites/sans.

Visit the repository at https://github.com/neutrons/drtsans.

Usage from provided front-ends
------------------------------

For end users go to
`QA version <http://scse-ui.ornl.gov:8080/>`_

Use `jupyter <https://jupyter.sns.gov/>`_ to have a play with the code.
The kernel to select is ``sans at ...``.

One can run scripts directly on `analysis <https://analysis.sns.gov/>`_ cluster.
To do that, open a terminal and activate the desired conda environment. The options are:

* ``sans`` the latest stable release
* ``sans-qa`` the future stable release (to be tested right before the next iteration)
* ``sans-dev`` the latest development version

The easiest way to start an interactive ipython session is by running

.. code-block:: shell

   $ drtsans

adding ``--qa`` or ``--dev`` will start the qa or development version respectively.
The ``drtsans`` wrapper script launches ipython with the selected conda environment located in ``/opt/anaconda/envs/`
and deactivates the conda environment when the session ends.


One must have an XCAMS account to use either the jupyter kernel provided above.

-----------------------------------------------
Set-up for development in a virtual environment
