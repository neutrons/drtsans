.. _reduction_output:

Reduction Output
================



canSAS format
-------------

Collective Action for Nomadic Small-Angle Scatterers (canSAS) provides standards and tools for the small-angle scattering user community. See more at https://www.cansas.org

canSAS defines a data format used in `drtsans`. During biosans, gpsans, and eqsans workflows, a canSAS file is generated that describes the 2-dimensional workspace I(Qx, Qy) with error.

Generating reports
------------------
Customized reports could be generated form reduction hdf log files. The script ``generate_reports`` extracts parameters specified in the generate_report.yaml file and prints them.

.. code-block:: bash

   generate_report  /path/to/hdf/log/file.hdf [/path/to/my/generate_report.yaml]

The last parameter here is a path to a users' YAML file with report parameters. If not provided, the [default parameters](https://github.com/neutrons/drtsans/blob/next/scripts/generate_report.yaml) will be used to create the report.

The yaml file contains the keys to extract from the hdf log file and short aliases to be used in the report

.. code-block:: bash

    options:
        reduction_information/drtsans/version : "drtsan version"
        reduction_information/mantid/version : "mantid version"
        reduction_information/reduction_parameters/background/runNumber : "background run number"
        reduction_information/special_parameters/background_transmission/value : "background_transmission"
        reduction_information/special_parameters/background_transmission/error : "background_transmission error"
        reduction_information/reduction_parameters/beamCenter/runNumber: "beam center run number"
        reduction_information/special_parameters/beam_center/x : "beam center x"
        reduction_information/special_parameters/beam_center/y : "beam center y"
        reduction_information/special_parameters/sample_transmission/value : "sample transmission"
        reduction_information/special_parameters/sample_transmission/error : "sample transmission error"

This yaml file generates the following output:

.. code-block:: bash

        drtsan version                                                   1.10.2+d20231019
        mantid version                                                              6.8.0
        background run number                                                       90335
        background_transmission                                        0.9254989091541436
        background_transmission error                               0.0029664327599940666
        beam center run number                                                      90277
        beam center x                                               0.0062761730094059885
        beam center y                                               -0.019878669581987165
        sample transmission                                            0.8939763002876431
        sample transmission error                                    0.002889993726769822

GRASP
-----
The command ``grasp_cg2`` can be used to save the pixel intensities from a gpsans reduction run into a GRASP-formatted ASCII file which can be loaded with the GRASP software package. Must be connected to analysis.sns.gov.

.. code-block:: bash

    >>> conda activate drtsans
    >>> grasp_cg2 -h

    usage: grasp_cg2 [-h] datafile output_directory

    drtsans GRASP. Must be connected to analysis.sns.gov

    positional arguments:
        datafile          The datafile to be processed.
        output_directory  The output directory that the GRASP file will be saved to.

    options:
        -h, --help        show this help message and exit
