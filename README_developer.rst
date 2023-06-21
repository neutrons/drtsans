.. _my-reference-label: Readmedev

========
drt-sans
========

Data Reduction Toolkit for Small Angle Neutron Scattering

This packages is a collection of functionality for reducing SANS data developed in collaboration with the instrument scientists at the High Flux Isotope Reactor (HFIR) and Spallation Neutron Source (SNS) at Oak Ridge National Laboratory.
While much of the functionality is generic, this implementation is aimed at reducing data from BIOSANS, EQSANS, and GPSANS.
As appropriate, this work is an abstraction layer on top of the mantid project.

**This is a python3 only package.**

-----------------------------------------------
Set-up for development in a virtual environment
-----------------------------------------------

These are the instructions for a person who is developing *both*
drt-sans and mantid in tandem. It assumes that you are using a virtual
environment and have a local build of mantid. As one caveat, for this
method, mantid must be build against the same version of python being
used in the virtual environment

1. Checkout the code.

  .. code-block:: shell

     $ git clone git@code.ornl.gov:sns-hfir-scse/sans/sans-backend.git
     $ cd sans-backend


2. Create the development environment.

Create a conda environment and activate it.

* It will use Mantid’s conda package.

* miniconda/conda is required.

* Linux command: Use the `libmamba solver <https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community>`_ for faster install

 .. code-block:: shell

    $ conda config --add channels conda-forge
    $ conda config --add channels mantid
    $ conda env create # will use environment.yml
    $ conda activate drtsans-dev

 **Warning**: if you create the conda environment on SNS’s analysis cluster, avoid naming your environment as ‘drtsans’,
 ‘drtsans-qa’ and ‘drtsans-dev’, which are reserved.

To make it easier to use, configure `direnv <https://direnv.net>`_ to activate the conda environment when you enter the directory.

3. Install the code in ``develop`` mode.

 .. code-block:: shell

    $ python setup.py develop

4. Try it out. Start ``python`` and try

 .. code-block:: python

    import mantid
    import drtsans

 Verify you can run the unit tests:

 .. code-block:: shell

    $ python -m pytest tests/unit/new/
    $ python -m pytest tests/integration/new/

 Some unit and integration tests require testing data.  Thus these tests can be either run on
 * SNS analysis cluster or
 * Local computer with mount to `/SNS/` archive. Here is the `instruction <https://code.ornl.gov/pf9/sns-mounts>`_ to mount SNS data directories.

5. When done, deactivate the virtual environment using

 .. code-block:: shell

    $ deactivate

 for virtual environment.  Or

 .. code-block:: shell

    $ deactivate

 for conda environment.


-----------------
Running the tests
-----------------
.. _running_tests:

The tests for this project are all written using `pytest <https://docs.pytest.org/en/latest>`_.
The `build pipeline <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/.gitlab-ci.yml>`_ currently `runs the unit tests and integration tests <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/test_job.sh>`_ separately using

.. code-block:: shell

   $ python -m pytest tests/unit/new/
   $ python -m pytest tests/integration/new/

This is one of the ways `pytest allows for selecting tests <https://docs.pytest.org/en/latest/usage.html#specifying-tests-selecting-tests>`_.
Specifying a directory or file will run all tests within that directory (recursively) or file.
Specifying a regular expression using ``-k`` will select all tests that match the regular expression independent of where they are defined

.. code-block:: shell

   $ python -m pytest -k test_samplelogs

To run an individual test within an individual file add ``::`` to the filename to specify the test

.. code-block:: shell

   $ python -m pytest tests/unit/new/drtsans/tof/eqsans/test_beam_finder.py::test_center_detector


--------------------------
Building the documentation
--------------------------
.. _building_docs:

The site can be built directly using

.. code-block:: shell

   $ sphinx-build -b html docs/ build/sphinx/html

or

.. code-block:: shell

   $ python setup.py build_sphinx

--------------------------
Installing the pre-commit hook
--------------------------

To automatically run the `pre-commit <https://pre-commit.com>`_ steps (e.g. linting) when adding a commit, install the git pre-commit hook.

.. code-block:: shell

   $ pre-commit install

To run pre-commit locally without committing

.. code-block:: shell

   $ pre-commit run --all

------------
Contributing
------------

Contributing is done through merge requests of code that you have the permission to add.
See `CONTRIBUTING.rst <CONTRIBUTING.rst>`_ for more information.

-----------------------------
Test Driven Development (TDD)
-----------------------------


* Test driven Development

   drtSANS development follows `test-driven development <https://en.wikipedia.org/wiki/Test-driven_development>`_ (TDD) process [1].
   All software requirements for SANS data reduction shall be converted to test cases before software is fully developed.
   All software developments are tracked by repeatedly testing the software against all test cases.

* Unit test

  All methods and modules shall have unit tests implemented.
  Unit tests are located in `repo/tests/unit/new <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/tree/next/tests/unit/new>`_.
  A unit test shall be created in the corresponding directory to the method or module that it tests against.

  Examples:

  * `drtsans/resolution.py <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/blob/next/drtsans/resolution.py>`_ and `tests/unit/new/drtsans/test_resolution.py <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/blob/next/tests/unit/new/drtsans/test_resolution.py>`_.
  * `drtsans/tof/eqsans/incoherence_correction_q1d.py <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/blob/next/drtsans/tof/eqsans/incoherence_correction_1d.py>`_ and `tests/unit/new/drtsans/tof/eqsans/test_incoherence_correction_q1d.py <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/blob/next/tests/unit/new/drtsans/tof/eqsans/test_incoherence_correction_q1d.py>`_.

* Integration test

  Integration test will test the combination of Individual modules and methods.
  Integration tests can be

  * general for all instrument, for instance `tests/integration/new/drtsans/test_stitch.py`.
  * specific to a suite of similar instruments, for instance `tests/integration/new/drtsans/mono/test_transmission.py` for all mono-wavelength instruments including Bio-SANS and GP-SANS.
  * specific to an individual instrument, for instance, `tests/integration/new/drtsans/mono/gpsans/test_find_beam_center.py` for GP-SANS and
    `tests/integration/new/drtsans/tof/eqsans/test_apply_solid_angle.py` for EQ-SANS.

* Testing data location

  Testing data are located on SNS data archive: `/SNS/EQSANS/shared/sans-backend/data`.

  Testing data for specific instruments have specific locations:

  - EQSANS: `/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/sns/eqsans/`
  - Bio-SANS: `/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/biosans/`
  - GP-SANS: `/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/gpsans/`

  Data files are referenced in the tests via the reference_dir pytest fixture.
  For instance, reference_dir.new.eqsans points to /SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/sns/eqsans/


------------------
Required libraries
------------------

* numpy: https://numpy.org/

* Mantid: https://www.mantidproject.org/, https://github.com/mantidproject/mantid

* Others: h5py, docutils, jsonschema, lmfit, matplotlib, mpld3, numexpr, pandas, sortedcontainers, tinydb, ipywidgets

* For unit and integration tests: pytest, pytest-xdist

* For documentation: sphinx, sphinxcontrib-napoleon,

* For linting and formatting: ruff which is configred in pre-commit
