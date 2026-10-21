.. _user.sensitivity_biosans:

BIOSANS sensitivity
===================

BIOSANS has a main detector and curved detector panels. The BIOSANS preparer
supports the main detector, wing detector, and midrange detector components;
the component is selected when constructing the preparer. The existing
template is ``scripts/prepare_sensitivities_biosans.py``.

Detector components
-------------------

Select one component for each sensitivity map:

``detector1``
   The main detector response normalization.

``wing_detector``
   The curved wing detector response normalization.

``midrange_detector``
   The optional curved midrange detector response normalization when present in
   the input data.

Curved panels are handled separately from the main detector. The preparer
detects the available curved panels from the first flood run and masks panels
that are not part of the selected component. This prevents counts from a
different panel from entering the selected component's sensitivity map.

Inputs and options
------------------

The common flood, direct-beam, dark-current, mask, beam-center-radius,
pixel-calibration, component-scaling, solid-angle, threshold, and preparation
method options are available. BIOSANS additionally supports a main-detector
mask angle through ``set_masks``. This is useful when a beam-trap or angular
mask is more appropriate than a circular beam-center mask.

BIOSANS also supports transmission correction:

* ``set_transmission_correction`` supplies the transmission flood and reference
  runs and can set the beam-trap size factor; and
* ``set_theta_dependent_correction_flag`` selects whether the transmission is
  applied as a function of detector angle.

Transmission runs should be paired with the flood run for which they were
collected. The preparer applies the configured masks and detector-component
handling while loading these auxiliary data.

Preparation behavior
--------------------

The direct-beam workspace is used to determine the beam center. The sensitivity
normalization is calculated from the flood workspace, not from the beam-center
workspace. For the selected component, the flood response is normalized,
thresholded, and exported as a processed NeXus sensitivity workspace.

Choose the moving-detector method when the flood measurements provide detector
coverage through motion. Otherwise use the detector-patch method and select
thresholds appropriate for the measured response distribution.

Component example
-----------------

.. code-block:: python

   from drtsans.mono.biosans.prepare_sensitivities_correction import (
       PrepareSensitivityCorrection,
   )

   preparer = PrepareSensitivityCorrection(component="wing_detector")
   preparer.set_flood_runs((12345,))
   preparer.set_direct_beam_runs((12346,))
   preparer.set_masks(None, "1-8,249-256", main_det_mask_angle=0.5)
   preparer.set_beam_center_radius(140.0)
   preparer.set_solid_angle_correction_flag(True)
   preparer.execute(False, 0.5, 2.0, "biosans_wing_sensitivity.nxs")

Use the production template for transmission, dark-current, calibration, and
component-specific settings.
