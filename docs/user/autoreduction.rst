.. _user.autoreduction:


Autoreduction
=============

Python scripts for automatic reduction of SANS data are provided in the ``scripts/autoreduction/`` directory.
These scripts are designed to be run by the SNS web monitor but also from the command line.

Required arguments for the scripts include:

- Events Nexus file, e.g., ``/SNS/EQSANS/IPTS-20196/nexus/EQSANS_89157.nxs.h5``
- Output directory, e.g., ``/SNS/EQSANS/IPTS-20196/shared/autoreduction/``


Example of the autoreduction service running the script:

.. code-block:: bash

    python reduce_EQSANS.py \
        /SNS/EQSANS/IPTS-20196/nexus/EQSANS_89157.nxs.h5 \
        /SNS/EQSANS/IPTS-20196/shared/autoreduction/

User can also run the script manually, for instance:

.. code-block:: bash

    python reduce_EQSANS.py \
        /SNS/EQSANS/IPTS-20196/nexus/EQSANS_89157.nxs.h5 \
        /tmp --no_publish

This call will not publish the results to the live data server,
but will save report file ``EQSANS_89157.html`` in the ``/tmp`` directory.
User can point the web browser to the generated report file file:///tmp/EQSANS_89157.html.

Input Configuration
-------------------
The autoreduction script will automatically look for an input configuration file in three locations,
in the following order of precedence:

1. Under the ``shared`` directory of the current IPTS and taking into account the run number.
   For example, for IPTS-20196 and run number 89157,
   the script will look for ``/SNS/EQSANS/IPTS-20196/shared/autoreduce/89157/reduction_options_89157.json``.
2. Under the ``shared`` directory of the current IPTS. For example, for IPTS-20196,
   the script will look for ``/SNS/EQSANS/IPTS-20196/shared/autoreduce/reduction_options.json``.
3. File ``/SNS/EQSANS/shared/autoreduce/reduction_options.json``.

Output Files
------------
For IPTS-20196 and run number 89157, all reduction files are saved under directory
``/SNS/EQSANS/IPTS-20196/shared/autoreduction/89157/``. Among them:

1. ``EQSANS_89157.html`` - the same report that shows up in the web monitor.
2. ``EQSANS_89157_Iq.dat`` and ``EQSANS_89157_Iq.png`` - 1D reduced data file and plot.
3. ``EQSANS_89157_Iqxqy.dat``, ``EQSANS_89157_Iqxqy.h5``, ``EQSANS_89157_Iqxqy.png`` - 2D reduced data files and plot.
4. ``EQSANS_89157_processed.nxs`` - I(lambda) per detector pixel.
5. ``EQSANS_162568_reduction_log.hdf`` - HDF5 file with reduction-derived parameters (transmission, center, etc).
6. ``autoreduce_89157.log`` - log file with detailed information about the reduction steps at level ``INFO`` and above.
7. ``reduction_options_89157.json`` - input configuration options, completed with any missing default values.
   If User requests to autoreduce the same run again, this file will be used as input configuration.

Resubmitting Autoreduction
--------------------------
User can resubmit the autoreduction for a given run,
perhaps after modifying one of the reduction options files under the current IPTS,
then clicking the `reduction` link located at the bottom on the web monitor page for that run.
For example:

.. figure:: /user/media/reduction_button.png
   :alt: click the reduction link
   :width: 500px


TroubleShooting
---------------
If autoreduction fails, any log error messages as well as the traceback will be saved in the HTML report file
as well as published to the live data server (if option ``--no_publish`` is not used).


Livereduction
=============

Python scripts for automatic reduction of SANS data are provided in the ``scripts/livereduction/`` directory.
These scripts are meant to be run automatically by the live reduction server, never directly by the User.
