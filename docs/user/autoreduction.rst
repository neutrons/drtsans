.. _user.autoreduction:


Autoreduction
=============

Python scripts for automatic reduction of SANS data are provided in the `scripts/autoreduction` directory.
These scripts are designed to be run from the command line and can be used to process data without manual intervention.

Required arguments for the scripts include:

- Events Nexus file, e.g., `/SNS/EQSANS/IPTS-20196/nexus/EQSANS_89157.nxs.h5`
- Output directory, e.g., `/SNS/EQSANS/IPTS-20196/shared/autoreduction/`


Example of the autoreduction service running the script:

code-block:: bash

    python reduce_EQSANS.py \
        /SNS/EQSANS/IPTS-20196/nexus/EQSANS_89157.nxs.h5 \
        /SNS/EQSANS/IPTS-20196/shared/autoreduction/

User can also run the script manually, for instance:

code-block:: bash

    python reduce_EQSANS.py \
        /SNS/EQSANS/IPTS-20196/nexus/EQSANS_89157.nxs.h5 \
        /tmp --no_publish --report_file EQSANS_89157.html

This call will generate the report file `EQSANS_89157.html` in the `/tmp` directory,
and will not publish the results to the live data server.
User can point the web browser to the generated report file file:///tmp/EQSANS_89157.html.
