=========================
:mod:`drtsans.resolution`
=========================

.. contents::

Smearing Pixels
===============

There are three different sources for smearing pixel width (X) and height (Y), ranked by their priority:

1. Reduction parameters `smearingPixelSizeX` and `smearingPixelSizeY`,
2. Barscan and tube-width pixel calibrations, and
3. Instrument definition file

A variety of scenarios giving rise to different **final smearing pixel sizes** are presented:

1. If no reduction parameters and no pixel calibration is supplied, then the instrument definition file
provides a smearing width :math:`w_0` and height :math:`h_0` for all pixels.

2. If reduction parameters `smearingPixelSizeX` and `smearingPixelSizeY` are supplied, and no
pixel calibration is supplied, then `smearingPixelSizeX` and `smearingPixelSizeY` become
the smearing width and height for all pixels.

3. If no reduction parameters are supplied but pixel calibration is supplied, then the
smearing width and height are taken from pixel calibration pixel sizes.

4. Finally, if reduction parameters `smearingPixelSizeX` and `smearingPixelSizeY` are supplied,
and a pixel calibration is also supplied, the smearing width :math:`w_i` of pixel :math:`i`
becomes

.. math::

    w_i = \frac{smearingPixelSizeX}{w_0} \cdot w_{pc, i},

where :math:`w_{pc, i}` is the pixel width of pixel :math:`i` provided by the pixel calibration.
An analogous relation follows for the final smearing height.

Here is a **diagram** of the functions involved in porting input reduction parameters `smearingPixelSizeX`
and `smearingPixelSizeY` into the function calculating the undeterminacy in momentum transfer:

.. graphviz::

   digraph foo {
      A1 [label="drtsans.mono.biosans.api.load_all_files", shape=box, href="#drtsans.mono.biosans.api.load_all_files"]
      A2 [label="drtsans.mono.gpsans.api.load_all_files", shape=box, href="#drtsans.mono.gpsans.api.load_all_files"]
      A3 [label="drtsans.tof.eqsans.api.load_all_files", shape=box, href="#drtsans.tof.eqsans.api.load_all_files"]
      B1 [label="drtsans.mono.meta_data.set_meta_data", shape=box, href="#drtsans.mono.meta_data.set_meta_data"]
      B2 [label="drtsans.tof.eqsans.meta_data.set_meta_data", shape=box, href="#drtsans.tof.eqsans.meta_data.set_meta_data"]
      C [label="drtsans.geometry.logged_smearing_pixel_size", shape=box, href="#drtsans.geometry.logged_smearing_pixel_size"]
      D1 [label="drtsans.mono.momentum_transfer.retrieve_instrument_setup", shape=box, href="#drtsans.mono.momentum_transfer.retrieve_instrument_setup"]
      D2 [label="drtsans.tof.eqsans.momentum_transfer.retrieve_instrument_setup", shape=box, href="#drtsans.tof.eqsans.momentum_transfer.retrieve_instrument_setup"]
      E [label="drtsans.resolution.InstrumentSetUpParameters", shape=box, href="#drtsans.resolution.InstrumentSetupParameters"]
      F [label="drtsans.resolution.calculate_sigma_theta_geometry", shape=box, fontcolor=blue, href="#drtsans.resolution.calculate_sigma_theta_geometry"]
      A1 -> B1
      A2 -> B1
      A3 -> B2
      B1 -> C
      B2 -> C
      C -> D1 -> E
      C -> D2 -> E
      E -> F;
   }


API
===

.. automodule:: drtsans.resolution
   :members:
   :private-members:
   :special-members:

