.. _user.sensitivity_eqsans:

EQSANS sensitivity
==================

The user-facing entry point is the ``prepare_sensitivities_eqsans`` command.
It reads an experiment-specific TOML configuration file and writes a processed
NeXus sensitivity file. The output path is a required command-line argument
and is deliberately not part of the configuration file.

In a development checkout, run the source script from the repository root:

.. code-block:: shell

   pixi run python scripts/prepare_sensitivities_eqsans.py \
      eqsans_sensitivity.toml \
      /path/to/sensitivity.nxs

When the conda package is installed, as in environments `sans_dev`, `sans_qa`, and `sans`,
invoke the packaged command directly:

.. code-block:: shell

   prepare_sensitivities_eqsans \
      eqsans_sensitivity.toml \
      /path/to/sensitivity.nxs

The configuration file uses uppercase keys matching the settings described
below. Settings omitted from the file retain their defaults. TOML does not
have a ``None`` value, so omit optional settings when the corresponding
correction is not needed. Invalid values and unknown keys are rejected.

For example:

.. code-block:: toml

   FLOOD_RUNS = [188982]
   DIRECT_BEAM_RUNS = [188969]

   BLOCKED_BEAM_RUNS = 188984
   EMPTY_CELL_RUNS = 188978

   TRANSMISSION_FLOOD_RUNS = 188974
   TRANSMISSION_EMPTY_CELL_RUNS = 188970
   TRANSMISSION_REFERENCE_RUNS = 188969

   SAMPLE_OFFSET = 285.0
   DETECTOR_OFFSET = 66.714
   LOW_TOF_CLIP = 1650.0
   HIGH_TOF_CLIP = 3150.0

   MASKED_PIXELS = "1-15,242-256"
   MASK_BEAM_CENTER_RADIUS = 40.0
   SOLID_ANGLE_CORRECTION = true
   SCALE_COMPONENTS = { detector1 = [1.004124, 1.057996, 1.0] }

   MIN_THRESHOLD = 0.1
   MAX_THRESHOLD = 3.0

The command prepares the ``detector1`` component. The settings are grouped
below in the same order as they appear in the source script.

Configuration settings
----------------------

Each run setting contains exactly one run: one integer, one string, or a
one-element TOML array containing either form. Numeric strings are converted to integers. Other strings are preserved as
run identifiers or NeXus paths. For example, all of the following are valid:
``12345``, ``"12345"``, ``[12345]``, and ``["12345"]``.

``FLOOD_RUNS``
   Flood measurement used to calculate the detector response. Required.

``DIRECT_BEAM_RUNS``
   Measurement used to determine the beam center. The beam center is used to
   center the detector and to mask the beam-stop region. Required.

``DARK_CURRENT_RUNS``
   Optional dark-current measurement. Omit this setting when no dark-current
   subtraction is required. Default: omitted.

``BLOCKED_BEAM_RUNS``
   Optional blocked-beam measurement. Omit this setting when this background
   is not required. Default: omitted.

``EMPTY_CELL_RUNS``
   Optional empty-cell measurement. Omit this setting when this background is
   not required. Default: omitted.

The EQSANS sensitivity is computed with the detector-patch method, which uses a
single flood image. The command therefore accepts exactly one flood run, and
rejects any run setting that contains more than one run.

Transmission settings
---------------------

Configure all three settings together to apply wavelength-dependent
transmission. Omitting all three disables this correction:

``TRANSMISSION_FLOOD_RUNS``
   Flood measurement used to determine transmission.

``TRANSMISSION_EMPTY_CELL_RUNS``
   Empty-cell transmission measurement.

``TRANSMISSION_REFERENCE_RUNS``
   Empty-beam reference measurement.

Each setting contains exactly one run. Providing only some of the three
settings is invalid. Transmission correction also requires
``EMPTY_CELL_RUNS``, because the transmission is applied to both the flood and
the empty cell before the empty cell is subtracted.

Geometry and data preparation
-----------------------------

``SAMPLE_OFFSET`` and ``DETECTOR_OFFSET``
   Float offsets in millimeters applied to flood and auxiliary data loads.
   These values should match the EQSANS reduction configuration; they are not
   inferred from the Nexus logs by this command. Defaults: ``0.0``.

``LOW_TOF_CLIP`` and ``HIGH_TOF_CLIP``
   Non-negative time-of-flight clipping widths in microseconds, as in the
   EQSANS reduction configuration. ``LOW_TOF_CLIP`` trims the fast-neutron end
   of each transmitted wavelength band, measured from the band's minimum TOF;
   ``HIGH_TOF_CLIP`` trims the slow-neutron end, measured back from the band's
   maximum TOF. The two widths are independent, so ``LOW_TOF_CLIP`` may exceed
   ``HIGH_TOF_CLIP``. Defaults: ``500.0`` and ``2000.0``.

``UNIVERSAL_MASK``
   A mask file path or a TOML array of detector IDs. Omit this setting when no
   universal mask is required. A Mantid workspace cannot be supplied through
   TOML. Default: omitted.

``MASKED_PIXELS``
   A pixel-range string such as ``"1-18,239-256"``. Omit this setting when no
   additional pixels should be masked. Default: omitted.

``MASK_BEAM_CENTER_RADIUS``
   Beam-center mask radius as a float in millimeters. Default: ``65.0``.

``PIXEL_CALIBRATION``
   ``False`` to disable pixel calibration, ``True`` to use the default
   calibration, or a calibration-file path. Default: ``false``.

``SCALE_COMPONENTS``
   A TOML inline dictionary mapping component names to positive ``[X, Y, Z]``
   scale factors. For EQSANS, the supported component is ``detector1``. Omit
   this setting when no component scaling is required. This corresponds to
   ``configuration.scaleComponents`` in the reduction configuration. Default:
   omitted.

``SOLID_ANGLE_CORRECTION``
   Boolean selecting whether flood data are divided by pixel solid angle before
   the sensitivity map is formed. The default is ``true`` so the map
   stores relative detector response without the solid-angle dependence.

Sensitivity calculation and output
----------------------------------

``MOVING_DETECTORS``
   Must be ``false``. EQSANS sensitivities are calculated with the
   detector-patch method, which fills the beam-stop region of a single flood by
   fitting the surrounding pixels of each tube. The moving-detector method,
   which combines floods with the beam-stop shadow on different pixels, is
   available only for GPSANS and BIOSANS. The command rejects ``true``.
   Default: ``false``.

``MIN_THRESHOLD`` and ``MAX_THRESHOLD``
   Floating-point lower and upper limits for normalized detector response.
   Pixels outside this range are masked in the output map. Defaults: ``0.5``
   and ``2.0``.

Processing order
----------------

Flood, blocked-beam, empty-cell, and transmission workspaces are normalized by
proton charge. The beam-center run is not normalized because it is used only to
locate the beam center. No beam-flux spectrum is used for sensitivity
preparation. Flood, blocked-beam, empty-cell, and
transmission workspaces use the same configured geometry and TOF clipping widths. The
blocked-beam and empty-cell backgrounds are subtracted before wavelength bins
are summed.

When transmission correction is configured, the transmission flood and empty
cell are compared with the empty-beam reference for each wavelength. The
resulting ``T(lambda, theta)`` is applied to the flood and empty-cell data
before their background subtraction and wavelength integration. Transmission
loads do not apply solid-angle correction because that correction belongs to
the sensitivity flood geometry.

The EQSANS workflow deliberately differs from the ordinary sample-reduction
workflow in one important respect: a sample reduction can apply sensitivity
normalization to the transmission runs before calculating transmission,
whereas sensitivity preparation is producing the sensitivity map and cannot
apply that output map to its own transmission inputs. The map's provenance logs
record the transmission configuration and this known distinction: when
transmission correction is configured, the ``transmission_sensitivity_corrected``
log is ``0``.

Beam-center data
----------------

The direct-beam workspace is loaded with the beam-center preparation options
needed to determine the center. EQSANS sample/detector offsets and TOF cuts are
not applied to that beam-center load as flood corrections. The beam-center
workspace is therefore not background-subtracted, transmission-corrected,
wavelength-integrated as a flood, or used as the output sensitivity map. The
beam-center region is masked on the prepared flood workspace before the
sensitivity response is calculated.

Uncertainties and output
------------------------

The initial EQSANS count uncertainty follows the Poisson model,
``sigma = sqrt(N)``. Mantid propagates uncertainty through charge
normalization, background subtraction, geometry operations, and transmission
correction. Detector-patch fitting uses inverse measured uncertainty as its
weight and propagates the full covariance matrix of the fitted parameters to
the response uncertainty. This full-covariance treatment is intentional and
is a documented deviation from the independent-parameter shortcut described
in the original requirements.

The exported processed NeXus file retains nominal output geometry and records
the processing choices in sample logs, including the drtsans version, flood
and auxiliary runs, beam centers, offsets, TOF clipping widths, masking,
normalization, transmission, thresholds, and preparation method.
