=======================
Developer documentation
=======================

See `solid_angle_correction.py <drtsans/solid_angle_correction.py>`_ for propely documented code

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

     $ git clone git@github.com:neutrons/drtsans.git
     $ cd drtsans


2. Create the development environment.

Pixi is used for managing environments, dependencies, packaging, and task execution.

* It will use Mantid’s conda package.

* Pixi is required. Install Pixi with the following command:

  .. code-block:: shell
   
     $ curl -fsSL https://pixi.sh/install.sh | bash
   

* Create the evironment by running

  .. code-block:: shell
   
     $ pixi install

* Activate the environment by running
  
  .. code-block:: shell
  
     $ pixi shell

3. Install the code in ``develop`` mode.

 .. code-block:: shell

    $ pip install -e .

4. Try it out. Start ``python`` and try

 .. code-block:: python

    import mantid
    import drtsans

 Verify you can run the unit tests:

 .. code-block:: shell

    $ python -m pytest tests/unit/
    $ python -m pytest tests/integration/

 Some unit and integration tests require testing data.  Thus these tests can be either run on
 * SNS analysis cluster or
 * Local computer with mount to `/SNS/` archive. Here is the `instruction <https://code.ornl.gov/pf9/sns-mounts>`_ to mount SNS data directories.

5. When done, deactivate the virtual environment using

 .. code-block:: shell

    $ exit


---------------------------
Git submodule for test data
---------------------------

To get the current version of the ``git-lfs`` associated with the currently checked out branch

.. code-block:: shell

   $ git submodule update --init

where ``update`` checks out the associated refspec and ``--init`` tells git to initialize the submodule if it isn't already.
To update the submodule to the latest commit on the branch being tracked

.. code-block:: shell

   $ git submodule update --remote --merge

where ``--remote`` tells git to fetch the latest changes from the upstream repository and ``--merge`` tells git to merge the changes into the working tree.

Then, to get the data files associated with the submodule, run

.. code-block:: shell

   $ cd tests/data/drtsans-data
   $ git lfs fetch

See the `git-submodules documentation <https://git-scm.com/book/en/v2/Git-Tools-Submodules>`_ for more detailed information.

-----------------
Running the tests
-----------------
.. _running_tests:

The tests for this project are all written using `pytest <https://docs.pytest.org/en/latest>`_.
The `build pipeline <https://github.com/neutrons/drtsans/blob/next/.github/workflows/test.yml>`_ currently `runs the unit tests and integration tests separately using

.. code-block:: shell

   $ python -m pytest tests/unit/
   $ python -m pytest tests/integration/

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

   $ pixi run build-docs

------------------------------
Installing the pre-commit hook
------------------------------

To automatically run the `pre-commit <https://pre-commit.com>`_ steps (e.g. linting) when adding a commit, install the git pre-commit hook.

.. code-block:: shell

   $ pre-commit install

To run pre-commit locally without committing

.. code-block:: shell

   $ pre-commit run --all


-----------------------------
Test Driven Development (TDD)
-----------------------------


* Test driven Development

   drtSANS development follows `test-driven development <https://en.wikipedia.org/wiki/Test-driven_development>`_ (TDD) process [1].
   All software requirements for SANS data reduction shall be converted to test cases before software is fully developed.
   All software developments are tracked by repeatedly testing the software against all test cases.

* Unit test

  All methods and modules shall have unit tests implemented.
  Unit tests are located in `repo/tests/unit/new <https://github.com/neutrons/drtsans/blob/next/tests/unit/>`_.
  A unit test shall be created in the corresponding directory to the method or module that it tests against.

  Examples:

  * `drtsans/resolution.py <https://github.com/neutrons/drtsans/blob/next/src/drtsans/resolution.py>`_ and `tests/unit/drtsans/test_resolution.py <https://github.com/neutrons/drtsans/blob/next/tests/unit/drtsans/test_resolution.py>`_.
  * `drtsans/tof/eqsans/incoherence_correction_1d.py <https://github.com/neutrons/drtsans/blob/next/src/drtsans/tof/eqsans/incoherence_correction_1d.py>`_ and `tests/unit/drtsans/tof/eqsans/test_incoherence_correction_1d.py <https://github.com/neutrons/drtsans/blob/next/tests/unit/drtsans/tof/eqsans/test_incoherence_correction_1d.py>`_.

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

--------------------------------------------
Package Build and Installation Instructions
--------------------------------------------

Python wheel

 .. code-block:: shell

    $ python -m build --wheel --no-isolation
    $ check-wheel-contents dist/drtsans-*.whl

Conda package

 .. code-block:: shell

    # create a conda package
    $ pixi run conda-build

.. _devdocs-standardnames:

------------------------------------------
Creating Candidate and Production Releases
------------------------------------------
- Follow the `Software Maturity Model <https://ornl-neutrons.atlassian.net/wiki/spaces/NDPD/pages/23363585/Software+Maturity+Model>`_ for continous versioning as well as creating release candidates and stable releases.
- Update file `release_notes.rst` with major fixes, updates and additions since last stable release.

--------------------------------------------
Standard variable names and what they denote
--------------------------------------------
* ``wavelength`` (was ``wl``) - wavelength
* ``sample_det_cent`` (was ``s2p``) - sample to detector center distance
* ``l1`` - source to sample distance
* **need definitions of** ``l2``, ``r1``, ``r2``, ``dwl``, ... These should likely get different names with unambigious meanings
* see e.g. solid_angle_correction.py for properly documented/named code


--------------------------------------
Random points about coding conventions
--------------------------------------

* All code will follow the `pep8 standard <https://www.python.org/dev/peps/pep-0008/>`_
* Docstrings will use `numpy formatting <https://numpydoc.readthedocs.io/en/latest/format.html>`_
* Code that uses Mantid algorithms will have an addional "Mantid algorithms used" section with links using inter-sphinx
* All internal functions should have an underbar, as is python standard
* Need more comments in meat of code explaining intent. For example, more clearly named variables (or leave the variables the same but give them better documentation)
* The development team still hasn't decided a standard for how errors/exceptions are handled
