.. _user.sensitivity_eqsans:

EQSANS sensitivity
==================

The user-facing entry point is the editable template
``scripts/prepare_sensitivities_eqsans.py``. Set the constants at the top of
that file for the experiment, then run the script in an environment configured
for Mantid and drtsans. The script constructs the preparer, applies the
configured settings, and writes the processed NeXus sensitivity file.

In a development checkout, run the template from the repository root:

.. code-block:: shell

   python scripts/prepare_sensitivities_eqsans.py

When the conda package is installed, it also provides the
``prepare_sensitivities_eqsans`` command. The command invokes the packaged
template:

.. code-block:: shell

   prepare_sensitivities_eqsans

The packaged template is installed under the drtsans package, while the source
template remains under the repository's ``scripts`` directory. To customize
the settings in an installed environment, copy the source template to a local
file, edit its settings, and run that file instead:

.. code-block:: shell

   python /path/to/prepare_sensitivities_eqsans.py

The copied template uses the installed drtsans package through its import and
does not need to be placed inside the package installation.

The template prepares the ``detector1`` component. The settings are grouped
below in the same order as they appear in the script.

Input runs
----------

Run settings accept one integer run number, one string run identifier or NeXus
path, or a sequence of either. For example: ``12345``, ``"EQSANS_12345"``,
``"/path/to/run.nxs"``, or ``(12345, 12346)``.

``FLOOD_RUNS``
   Required flood measurements used to calculate the detector response.

``DIRECT_BEAM_RUNS``
   Measurements used to determine the beam center. Provide one run or one run
   per flood run.

``DARK_CURRENT_RUNS``
   Optional dark-current measurements. Use ``None`` when no dark-current
   subtraction is required.

``BLOCKED_BEAM_RUNS``
   Optional blocked-beam measurements. Use ``None`` when this background is not
   required.

``EMPTY_CELL_RUNS``
   Optional empty-cell measurements. Use ``None`` when this background is not
   required.

An auxiliary run setting may contain one run, which is broadcast to all flood
runs, or one run per flood run. The script validates these pairings.

Transmission settings
---------------------

Configure all three settings together to apply wavelength-dependent
transmission:

``TRANSMISSION_FLOOD_RUNS``
   Flood measurements used to determine transmission.

``TRANSMISSION_EMPTY_CELL_RUNS``
   Empty-cell transmission measurements.

``TRANSMISSION_REFERENCE_RUNS``
   Empty-beam reference measurements.

Each setting accepts one run or one run per flood run. Leave all three as
``None`` to disable this transmission calculation.

Geometry and data preparation
-----------------------------

``SAMPLE_OFFSET`` and ``DETECTOR_OFFSET``
   Float offsets in millimeters applied to flood and auxiliary data loads.
   These values should match the EQSANS reduction configuration; they are not
   inferred from the Nexus logs by this template.

``LOW_TOF_CLIP`` and ``HIGH_TOF_CLIP``
   Float lower and upper time-of-flight limits in microseconds.

``UNIVERSAL_MASK``
   A mask file path, Mantid mask workspace, list of detector IDs, or ``None``.

``MASKED_PIXELS``
   A pixel-range string such as ``"1-18,239-256"``, or ``None``.

``MASK_BEAM_CENTER_RADIUS``
   Beam-center mask radius as a float in millimeters.

``PIXEL_CALIBRATION``
   ``False`` to disable pixel calibration, ``True`` to use the default
   calibration, or a calibration-file path.

``SCALE_COMPONENTS``
   A dictionary mapping component names to positive ``[X, Y, Z]`` scale
   factors, or ``None``. For EQSANS, the supported component is ``detector1``.
   This corresponds to ``configuration.scaleComponents`` in the reduction
   configuration.

``SOLID_ANGLE_CORRECTION``
   Boolean selecting whether flood data are divided by pixel solid angle before
   the sensitivity map is formed. The template default is ``True`` so the map
   stores relative detector response without the solid-angle dependence.

Sensitivity calculation and output
----------------------------------

``MOVING_DETECTORS``
   Boolean selecting the moving-detector method (``True``) or detector-patch
   method (``False``).

``MIN_THRESHOLD`` and ``MAX_THRESHOLD``
   Floating-point lower and upper limits for normalized detector response.
   Pixels outside this range are masked in the output map.

``OUTPUT_SENSITIVITY_FILE``
   Output path for the processed NeXus sensitivity file.

Processing order
----------------

Each EQSANS workspace is normalized by proton charge. No beam-flux spectrum is
used for sensitivity preparation. Flood, blocked-beam, empty-cell, and
transmission workspaces use the same configured geometry and TOF limits. The
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
record the transmission configuration and this known distinction.

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
and auxiliary runs, beam centers, offsets, TOF limits, masking, normalization,
transmission, thresholds, and preparation method.
