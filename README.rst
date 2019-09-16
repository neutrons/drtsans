drt-sans
========

For end users go to `next version <http://scse-sans-demo.ornl.gov/>`_
`QA version <http://scse-ui.ornl.gov:8080/>`_

Use `jupyter <https://jupyter.sns.gov/>`_ to have a play with the
code. The kernel to select is ``sans at ...``.


**This is a python3 only package.**

Using the Docker packaged environment
-------------------------------------

This the instructions for someone who wants to use the Docker container
created through the automated build pipeline to develop drt-sans, use
drt-sans to develop reduction scripts, or test existing drt-sans
functionality. The SNS analysis cluster does not have Docker installed
and Docker is required to follow these instructions.

1. (If not installed) Install Docker: https://docs.docker.com/install/
2. Download the latest sans-backend-run.sh script from the feature, release, or master branch for which you are testing: http://scse-mantid-demo.ornl.gov/sans-backend
3. Run the script with ``sudo bash sans-backend-run.sh -h`` to see the help menu.

Current options include:
* -i) launches a bash shell
* -u) forces an update of the application.
* -h) prints this message.

You must download the wrapper script from the above link as the build process modifies the copy in version control.

Set-up for development in a virtual environment
-----------------------------------------------

This is the instructions for a person who is developing *both*
drt-sans and mantid in tandem. It assumes that you are using a virtual
environment and have a local build of mantid. As one caveat, for this
method, mantid must be build against the same version of python being
used in the virtual environment

1. Checkout the code.

.. code-block:: shell

   $ git clone git@code.ornl.gov:sns-hfir-scse/sans/sans-backend.git
   $ cd sans-backend


2. Create the virtual environment and activate it. The
   ``--system-site-packages`` argument lets it use things installed on
   the system as a whole for compiling mantid

.. code-block:: shell

   $ virtualenv -p $(which python3) --system-site-packages .venv
   $ source .venv/bin/activate

3. Add the built version of mantid to the python path in the virtual
   environment

.. code-block:: shell

   $ python ~/build/mantid/bin/AddPythonPath.py

4. Install the requirements for running the code

.. code-block:: shell

   $ pip install -r requirements.txt -r requirements_dev.txt

5. Install the code in ``develop`` mode.

.. code-block:: shell

   $ python setup.py develop

6. Try it out. Start ``python`` and try

.. code-block:: python

   import mantid
   import ornl

Verify you can run the unit tests:

.. code-block:: shell

   $ python -m pytest tests/unit/new/

7. Be done. Deactivate the virtual environment using

.. code-block:: shell

   $ deactivate

As an alternative, you can use `direnv <https://direnv.net>`_ to do a
fair amount of the work. Unfortunately, it doesn't currently handle
``--system-site-packages`` correctly so you'll have to work around
it. Create the virtual environment using

.. code-block:: shell

   $ virtualenv -p $(which python3) --system-site-packages .direnv/python-$(python3 -c "import platform as p;print(p.python_version())")

Then you create a file ``.envrc`` in your source directory

.. code-block:: shell

   $ echo "layout python3" > .envrc

Finally, tell direnv that you want it to work in this directory

.. code-block:: shell

   $ direnv allow

Then follow steps 3-6 from above. After this, the virtual environment
with load when you enter the source tree, and unload when you leave.

Running the tests
-----------------

The tests for this project are all written using `pytest <https://docs.pytest.org/en/latest>`_.
The `build pipeline <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/.gitlab-ci.yml>`_ currently `runs the unit tests and integration tests <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/test_job.sh>`_ separately using

.. code-block:: shell

   $ python -m pytest tests/unit/new/
   $ python -m pytest tests/integration/new/

This is one of the ways `pytest allows for selecting tests <https://docs.pytest.org/en/latest/usage.html#specifying-tests-selecting-tests>`_.
Specifying a directory or file will run all tests within that directory (recursively) or file.
Specifying a regular expression using ``-k`` will select all tests that match the regular expression independent of where they are defined

.. code-block:: shell

   $ python -m pytest -k momentum_transfer

To run an individual test within an individual file add ``::`` to the filename to specify the test

.. code-block:: shell

   $ python -m pytest tests/integration/new/ornl/sans/sns/eqsans/test_momentum_transfer.py::test_api


Building the documentation
--------------------------

The site can be build directly using

.. code-block:: shell

   $ sphinx-build -b html docs/source/ build/sphinx/html

or

.. code-block:: shell

   $ python setup.py build_sphinx
