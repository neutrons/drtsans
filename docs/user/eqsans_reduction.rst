.. _user.eqsans_reduction:

EQSANS Reduction
================

This page collects reduction behavior that is specific to EQSANS.

Monochromatic Mode
-------------------

EQSANS can be operated in monochromatic mode, for runs
where the chopper system is phased to transmit a typically narrow wavelength band,
instead of the broad time-of-flight spectrum used in normal operation.
This mode is recorded in the raw data file through the process variable ``BL6:Chop:Skf16:MCON``,
set by the instrument's data acquisition system — it is not a setting in the reduction configuration.
This process variable has the associated alias ``MCON16``, which is the name under which it appears
in the sample logs of the raw data file, and thus the name `drtsans` looks up.
The behavior described below is triggered by the *value* of this boolean sample log; the mere
presence of the log in the file is not enough to trigger monochromatic mode.

When the ``MCON16`` sample log evaluates to ``True``, `drtsans` overrides the wavelength bin width
passed in the configuration file so that the entire transmitted band is treated as one bin spanning
the full range of wavelengths present in the data. Monochromatic mode is incompatible with
frame-skipping mode.

Effect on the time-of-flight clippings
++++++++++++++++++++++++++++++++++++++

Configuration options ``"cutTOFmin"`` and ``"cutTOFmax"`` discard neutrons at the two ends of
the time-of-flight frame, where the transmitted intensity tails off. Their default values, 500
and 2000 microseconds, are chosen for the broad wavelength band of normal operation. The
narrow band of monochromatic mode is comparable to, or smaller than, the span those clippings
remove, so applying them would discard most of the measured signal. For instance, run 177103
transmits the band 9.446-10.500 Angstrom; the default clippings amount to 0.140 Angstrom
below and 0.560 Angstrom above, leaving 9.586-9.940 Angstrom, barely a third of the band.
For a still narrower band the clippings exceed the band width altogether and the reduction
fails outright.

`drtsans` therefore ignores both clippings when the ``MCON16`` sample log evaluates to
``True``, keeping the whole transmitted band. This is not configurable: the values in the
reduction configuration file are left untouched but have no effect on a monochromatic run. The
override is reported in the reduction log:

.. code-block:: text

   Monochromatic mode detected: overriding TOF clippings (500.0, 2000.0) micro seconds with zero, to preserve the whole transmitted wavelength band

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

The clippings actually applied are recorded alongside it, under
``reduction_information/sample_logs/main/low_tof_clip`` and ``.../high_tof_clip``. Both read
zero for a monochromatic-mode reduction, whatever ``"cutTOFmin"`` and ``"cutTOFmax"`` were set
to in the reduction configuration.


Other EQSANS-specific documentation
------------------------------------

.. toctree::
   :maxdepth: 1

   /user/corrections/inelastic_incoherent
   /user/gpr_analysis
