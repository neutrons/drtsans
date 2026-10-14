.. _user.normalizations:

Normalizations
==============

Normalization scales measured detector data so that results can be compared
across runs, detector pixels, and instruments. The Master Requirements Document
describes normalization at several stages:

* pixel-area normalization, which accounts for the physical area of each
  detector pixel;
* detector-sensitivity normalization, which accounts for relative pixel
  response;
* exposure normalization by measurement duration, beam-monitor value or
  spectrum, or proton charge on target with a measured flux spectrum; and
* final absolute normalization by illuminated volume, an attenuated empty-beam
  measurement, or a calibrated standard sample.

These operations scale measured data; they do not represent a physical
background subtraction from the input data.

.. _user.sensitivity:

Sensitivity
-----------

Sensitivity is a detector-response normalization. A sensitivity map stores the
relative response of each detector pixel and the reduction workflow divides the
input data by that response.

When generating the sensitivity map, ``SOLID_ANGLE_CORRECTION`` determines
whether the input flood data are divided by the solid angle of each pixel before the
map is formed. When ``SOLID_ANGLE_CORRECTION = True``, the solid-angle
dependence is removed and the map stores only the relative pixel response.
When it is ``False``, the solid-angle dependence remains part of the map.

This option defaults to ``True`` in the preparation template and reduction
configuration because the reduction workflow applies solid-angle normalization
to the input data separately before dividing by the sensitivity map.

The common sensitivity-preparation workflow and the beamline-specific guides
are collected below:

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :caption: Beamline guides

   sensitivity_gpsans
   sensitivity_biosans
   sensitivity_eqsans

.. contents:: Sensitivity sections
   :local:
   :depth: 1

Common workflow
~~~~~~~~~~~~~~~

All three beamlines use the same high-level preparation sequence:

1. Load the flood and, when supplied, direct-beam measurements.
2. Apply the configured dark-current, mask, pixel-calibration, geometry, and
   solid-angle options.
3. Determine the beam center and mask the direct-beam region from the flood
   data. The beam-center workspace is used to determine the center; the flood
   workspace is the workspace from which the sensitivity normalization is
   calculated.
4. Normalize and combine the flood data.
5. Identify pixels outside the configured normalized-count thresholds and save
   the resulting map as a processed NeXus file.

The common implementation supports two preparation methods. The moving-detector
method combines measurements after mapping them to a common detector geometry.
The detector-patch method fits the response of each detector patch and retains
the detector geometry of the exported map. The available method depends on the
beamline data and should be selected to match the way the flood measurements
were acquired.

Uncertainty and provenance
~~~~~~~~~~~~~~~~~~~~~~~~~~

The preparation workflow keeps Mantid's propagated uncertainties with the
flood data. For event data, the initial count uncertainty follows the Poisson
model, ``sigma = sqrt(N)``, and subsequent subtraction, normalization,
geometry, and transmission operations propagate that uncertainty. The
detector-patch fit uses the measured uncertainties as weights and propagates
the full parameter covariance matrix into the fitted response. This is a
deliberate implementation choice for the detector-patch method rather than an
assumption that the fitted parameters are independent.

The exported sensitivity workspace contains processing logs. These include the
software version, preparation method, input runs, beam-center information,
masking and threshold options, geometry corrections, normalization method,
and beamline-specific corrections. These logs are intended to make a map
auditable and reproducible; they do not replace the original run metadata.

What sensitivity preparation does not do
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sensitivity preparation is not a sample reduction. It creates a response map
for later use and therefore does not divide the output sensitivity workspace by
itself. In particular, the beam-center workspace is used for beam-center
determination and is not treated as a sensitivity flood workspace.
Beamline-specific transmission treatment is described on the EQSANS and BIOSANS
pages because its available inputs and order of operations differ from the
common workflow.

Other normalization methods
----------------------------

The reduction workflows also normalize measurements by quantities such as
proton charge or monitor counts. The available method and its configuration are
beamline-specific. See the relevant reduction guide for those workflow options.
