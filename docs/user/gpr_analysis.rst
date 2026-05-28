.. _user.gpr_analysis:

Gaussian Process Regression Analysis
=====================================

The algorithm used here was developed by Changwoo Do at ORNL and has been optimized for neutron scattering data.

Automatic Use in Autoreduction
-------------------------------

GPR analysis runs automatically during EQSANS autoreduction. You don't need to do anything special.
The results appear in your autoreduction HTML report as interactive Plotly plots, and separate PNG and DAT files
are saved alongside your regular reduction outputs.

If you want to turn it off, add this to your reduction configuration:

.. code-block:: json

    {
        "configuration": {
            "enableGPR": false
        }
    }

Command-Line Usage
------------------

You can also run GPR analysis independently on any I(Q) data file:

.. code-block:: bash

    python -m drtsans.gpr input_file_Iq.dat

This will create two files in the same directory:

- ``input_file_Iq_gpr.png`` - Plot showing experimental data and GPR fit
- ``input_file_Iq_gpr.dat`` - ASCII file with Q, I_GPR, I_GPR_error, dQ

Specify a different output directory:

.. code-block:: bash

    python -m drtsans.gpr input_file_Iq.dat --output-dir /path/to/output

Interactive Terminal Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For an interactive experience, use the ``--tui`` flag to launch a terminal-based interface:

.. code-block:: bash

    python -m drtsans.gpr --tui

This presents a menu where you can select files and adjust parameters without typing long command lines.

Shell Tab Completion
~~~~~~~~~~~~~~~~~~~~

Enable tab completion for faster command-line use:

.. code-block:: bash

    # Add to your ~/.bashrc or ~/.zshrc:
    eval "$(register-python-argcomplete drtsans-gpr)"

Then you can tab-complete file paths and options when using the GPR command.

Python API
----------

For more control, use the Python API directly in your scripts or notebooks.

Basic Usage
~~~~~~~~~~~

.. code-block:: python

    from drtsans.gpr import generate_gpr_analysis
    from drtsans.dataobjects import IQmod

    # Load your I(Q) data
    iq_data = IQmod.read_csv("EQSANS_89157_Iq.dat")

    # Run GPR analysis
    html_report, png_files, dat_files = generate_gpr_analysis(
        iq_data,
        output_dir="/path/to/output",
        base_filename="EQSANS_89157"
    )

    print(f"Generated: {png_files[0]}")
    print(f"Data saved to: {dat_files[0]}")

Multiple Profiles (Wedges or Time Slices)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # If you have multiple I(Q) profiles from wedge analysis
    wedge_profiles = [iq_wedge_0, iq_wedge_1, iq_wedge_2, iq_wedge_3]

    html_report, png_files, dat_files = generate_gpr_analysis(
        wedge_profiles,
        output_dir="/path/to/output",
        base_filename="EQSANS_89157"
    )

    # Creates separate files for each wedge:
    # EQSANS_89157_wedge_0_Iq_gpr.png, EQSANS_89157_wedge_0_Iq_gpr.dat
    # EQSANS_89157_wedge_1_Iq_gpr.png, EQSANS_89157_wedge_1_Iq_gpr.dat
    # etc.

Custom Parameters
~~~~~~~~~~~~~~~~~

Adjust the GPR algorithm behavior with optional parameters:

.. code-block:: python

    html_report, png_files, dat_files = generate_gpr_analysis(
        iq_data,
        output_dir="/path/to/output",
        base_filename="EQSANS_89157",
        use_log_I=False,          # Log transform intensity? (default: False)
        use_log_Q=True,           # Log transform Q? (default: True)
        lmbda=0.25,               # Smoothness parameter (default: 0.25)
        background_filter_width=0.2,      # Background window (default: 0.2)
        background_intensity_offset=1.0   # Background offset (default: 1.0)
    )

File-Based Interface
~~~~~~~~~~~~~~~~~~~~

Process .dat files directly:

.. code-block:: python

    from drtsans.gpr import run_gpr_from_file

    png_path, dat_path = run_gpr_from_file(
        "EQSANS_89157_Iq.dat",
        output_dir="/path/to/output"
    )

Understanding the Parameters
-----------------------------

lmbda (Lambda)
~~~~~~~~~~~~~~

Controls the smoothness of the GPR fit. Smaller values give smoother curves:

- ``lmbda=0.1`` - Very smooth, might over-smooth real features
- ``lmbda=0.25`` - Default, good balance for most SANS data
- ``lmbda=0.5`` - Less smooth, follows data more closely
- ``lmbda=1.0`` - Minimal smoothing, close to the data

use_log_Q
~~~~~~~~~

Whether to apply log transformation to Q values:

- ``True`` (default) - Recommended for SANS data spanning multiple decades
- ``False`` - Use for linear Q-scale data

use_log_I
~~~~~~~~~

Whether to apply log transformation to intensity values:

- ``False`` (default) - Works well for most SANS profiles
- ``True`` - Can help with data spanning many orders of magnitude

background_filter_width
~~~~~~~~~~~~~~~~~~~~~~~

Sets the window size for iterative background estimation (as a fraction of total Q range):

- ``0.1`` - Narrow window, local background
- ``0.2`` (default) - Balanced background estimation
- ``0.5`` - Wide window, global background trends

Output Files
------------

GPR analysis creates two types of output files:

PNG Files
~~~~~~~~~

High-resolution plots (300 DPI) showing:

- Original experimental data points with error bars
- GPR fit curve (solid line)
- Uncertainty bands (shaded regions)
- Log-log axes for typical SANS visualization

These plots are suitable for publications and presentations.

DAT Files
~~~~~~~~~

ASCII text files with tab-separated columns:

.. code-block:: text

    # Q (1/A)    I_GPR (1/cm)    I_GPR_err (1/cm)    dQ (1/A)
    0.005000e+00    1.234567e+06    1.234567e+04    2.500000e-04
    0.005500e+00    1.123456e+06    1.123456e+04    2.750000e-04
    ...

These files can be loaded into other analysis software or plotting tools.

HTML Reports
~~~~~~~~~~~~

When running through autoreduction, GPR generates interactive Plotly plots embedded in the HTML report.
These allow you to:

- Zoom into regions of interest
- Hover over points to see exact values
- Toggle between experimental and fitted data
- Pan across the Q range

References
----------

**Algorithm**: Based on radial basis function (RBF) kernel GPR with log-space transformations

**Original Development**: Changwoo Do, Oak Ridge National Laboratory

**Integration**: May 2026 - Added to drtsans for EQSANS autoreduction

For more details on the statistical methods, see the API documentation: :ref:`drtsans.gpr module <api.drtsans.gpr>`
