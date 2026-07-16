.. _user.eqsans_reduction:

EQSANS Reduction
================

This page collects reduction behavior that is specific to EQSANS.

Monochromatic Mode
-------------------

EQSANS can be operated in monochromatic mode, for runs
where the chopper system is phased to transmit a typically narrow wavelength band,
instead of the broad time-of-flight spectrum used in normal operation.
This mode is recorded in the raw data file through the ``"monochromatic"`` sample log,
set by the instrument's data acquisition system — it is not a setting in the reduction configuration.
The behavior described below is triggered by the *value* of this boolean sample log; the mere
presence of the log in the file is not enough to trigger monochromatic mode.

When the ``"monochromatic"`` sample log evaluates to ``True``, `drtsans` overrides the wavelength bin width
passed in the configuration file so that the entire transmitted band is treated as one bin spanning
the full range of wavelengths present in the data. Monochromatic mode is incompatible with
frame-skipping mode.

Effect on wavelength-dependent corrections
++++++++++++++++++++++++++++++++++++++++++

The elastic reference normalization and inelastic incoherent compensation described in
:doc:`Wavelength-dependent corrections </user/corrections/inelastic_incoherent>` both work by
comparing the scattered intensity across several wavelength bins relative to a reference
wavelength. With a single wavelength bin, these corrections become meaningless.

To avoid this, `drtsans` automatically bypasses both corrections whenever a single wavelength
bin is detected, regardless of whether ``"fitInelasticIncoh"`` or ``"elasticReference"`` were
requested in the reduction configuration. Each bypassed correction is reported with a warning
in the reduction log:

.. code-block:: text

   Only 1 wavelength bin present in frame 0: bypassing elastic reference normalization correction.
   Only 1 wavelength bin present in frame 0: bypassing inelastic incoherence correction.

Verifying monochromatic-mode reduction
++++++++++++++++++++++++++++++++++++++

The wavelength bin width used for a given reduction is recorded in the reduction log HDF5
file, under ``reduction_information/sample_logs/main/wavelength_bin_width``. For a
monochromatic-mode reduction, this dataset contains a single value equal to the width of the
whole wavelength band, rather than one value per ``"wavelengthStep"``-sized bin, which can be
used to confirm that the monochromatic path was taken.


Other EQSANS-specific documentation
------------------------------------

.. toctree::
   :maxdepth: 1

   /user/corrections/inelastic_incoherent
   /user/gpr_analysis
