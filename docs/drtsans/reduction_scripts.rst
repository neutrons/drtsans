=================
Reduction Scripts
=================

.. contents::

The following python scripts can be used as the entry points for reduction of SANS data
for each instrument

- `scripts/biosans_reduction.py <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/blob/next/scripts/biosans_reduction.py>`_
- `scripts/eqsans_reduction.py <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/blob/next/scripts/eqsans_reduction.py>`_
- `scripts/gpsans_reduction.py <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/blob/next/scripts/gpsans_reduction.py>`_

These scripts receive as argument the path to a `*.json` file containing all necessary reduction parameters. In the
active `sans` conda environment, and assuming we are at the root of the drtsans repository:

.. code-block:: bash

   python ./scripts/biosans_reduction.py /path/to/my.json

Examples of reduction parameter files are available for each instrument, in the data repository

- `BIOSANS examples <https://code.ornl.gov/sns-hfir-scse/infrastructure/test-data/drtsans-data/-/blob/main/ornl/sans/hfir/biosans/reduction_parameters/README.md>`_
- `EQSANS examples <https://code.ornl.gov/sns-hfir-scse/infrastructure/test-data/drtsans-data/-/blob/main/ornl/sans/sns/eqsans/reduction_parameters/README.md>`_
- `GPSANS examples <https://code.ornl.gov/sns-hfir-scse/infrastructure/test-data/drtsans-data/-/blob/main/ornl/sans/hfir/gpsans/reduction_parameters/README.md>`_
