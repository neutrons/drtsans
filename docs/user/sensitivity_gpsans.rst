.. _user.sensitivity_gpsans:

GPSANS sensitivity
==================

GPSANS uses the common sensitivity preparer for the CG2 detector. The example
template is ``scripts/prepare_sensitivities_gpsans.py``.

Inputs and options
------------------

The usual preparation uses one or more flood runs and, when available, direct
beam runs. GPSANS flood data are normalized using the measured beam current.
The direct-beam data determine the beam center and the configured beam-center
region is removed from the flood data before the response is calculated.

The template exposes these common options:

* a universal mask and a detector-pixel list;
* a beam-center mask radius in millimeters;
* optional dark-current runs;
* optional pixel calibration;
* optional component scaling; and
* solid-angle correction.

Set ``MOVING_DETECTORS`` to ``True`` when the flood measurements cover the
detector by moving it. Set it to ``False`` to use the detector-patch method.
``MIN_THRESHOLD`` and ``MAX_THRESHOLD`` define the acceptable range of the
normalized response. Pixels outside that range are masked in the output map.

Preparation behavior
--------------------

GPSANS uses the common flood-processing path. The same configured mask,
calibration, geometry, and solid-angle choices are applied to the flood data
that contribute to the map. The beam-center data are loaded to locate the
beam, while the beam-center masking and sensitivity normalization are performed
on the flood workspace.

The output is a processed NeXus sensitivity workspace. The file name and
location are set by ``SENSITIVITY_FILE`` in the template. The reduction
workflow can then use this file to normalize sample data collected with the
same detector configuration.

Minimal example
---------------

The essential public API is:

.. code-block:: python

   from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection

   preparer = PrepareSensitivityCorrection(instrument="CG2", component="detector1")
   preparer.set_flood_runs((90091, 90093, 90095))
   preparer.set_direct_beam_runs((90092, 90094, 90096))
   preparer.set_masks(None, "1-12,239-256")
   preparer.set_beam_center_radius(140.0)
   preparer.set_pixel_calibration_flag(True)
   preparer.set_solid_angle_correction_flag(True)
   preparer.execute(True, 0.5, 2.0, "gpsans_sensitivity.nxs")

The template contains the complete set of editable constants and is preferred
when preparing a production map.
