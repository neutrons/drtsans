.. _user.eqsans_reduction:

EQSANS Reduction
================

This page collects reduction behavior that is specific to EQSANS.

The transmitted wavelength band
-------------------------------

The band of wavelengths reaching the sample is derived from the phases of the disk choppers. A
neutron of wavelength :math:`\lambda` leaves the moderator :math:`t_0(\lambda)` after the pulse
and reaches a chopper at distance :math:`D` at time :math:`t_0(\lambda) + D / v(\lambda)`; it is
transmitted when that time falls inside the aperture's opening window. Both edges of the band
therefore carry this emission-delay correction, which displaces the band towards shorter
wavelengths by an amount inversely proportional to the chopper's distance to the source,
typically 0.05 to 0.09 Angstrom.

The data acquisition system phases the choppers using the purely geometric calculation, which
ignores the emission delay, so the band that reaches the sample is displaced from the nominal
setpoint. Run 177103, for example, is phased for a round 9.500-10.500 Angstrom and transmits
9.446-10.410 Angstrom.

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
transmits the band 9.446-10.410 Angstrom; the default clippings amount to 0.140 Angstrom
below and 0.560 Angstrom above, leaving 9.586-9.849 Angstrom, barely a quarter of the band.
For a still narrower band the clippings exceed the band width altogether and the reduction
fails outright.

`drtsans` therefore ignores both clippings when the ``MCON16`` sample log evaluates to
``True``, keeping the whole transmitted band. This is not configurable: the values in the
reduction configuration file are left untouched but have no effect on a monochromatic run. The
override is reported in the reduction log:

.. code-block:: text

   Monochromatic mode detected: overriding TOF clippings (500.0, 2000.0) micro seconds with zero, to preserve the whole transmitted wavelength band

Cross-check against the requested band
++++++++++++++++++++++++++++++++++++++

The data acquisition system records the band it was asked to deliver in two more process
variables, ``MCWL16`` (the middle of the band, in Angstrom) and ``MCWLSpread16`` (the width of
the band, as a percent of the middle). The requested band therefore spans

.. math::

   \left[ \left(1 - \frac{p}{200}\right) \lambda_0, \left(1 + \frac{p}{200}\right) \lambda_0 \right]

for a center :math:`\lambda_0` and a spread :math:`p`. `drtsans` derives the band from the
chopper phases independently, and compares the two so that a mis-phased chopper set cannot pass
unnoticed. The figure of merit is the fraction of the *requested* band that the choppers are
phased for:

.. code-block:: text

   overlap = |requested band intersect geometric band| / |requested band|

The comparison uses the *geometric* band, computed from the chopper phases with no
emission-delay correction, because that is the quantity the data acquisition system phases the
choppers against. The band that actually reaches the sample is shifted towards shorter
wavelengths by the delayed emission of neutrons from the moderator, by up to 0.09 Angstrom, so
comparing the request against it would report a disagreement on every run. What the comparison
detects is chopper phases set for a band other than the one requested.

============================ ==========================================================
Overlap                      Behavior
============================ ==========================================================
90% or more                  reported at information level; the reduction proceeds
between 80% and 90%          reported as a warning; the reduction proceeds
less than 80%                reported as an error, and the reduction stops with a
                             ``ValueError``
============================ ==========================================================

Runs taken before these two process variables were introduced carry neither log, even when
flagged as monochromatic. Such a run cannot be checked; the omission is reported as a warning
and the reduction proceeds:

.. code-block:: text

   Monochromatic mode: cannot verify the transmitted wavelength band because sample log(s) MCWL16, MCWLSpread16 not found. Skipping the check.

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
