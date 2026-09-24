.. _user.corrections.attenuation:

Absolute Scaling with an Attenuated Direct Beam (GPSANS)
========================================================

This correction applies only to GPSANS. BIOSANS and EQSANS do not implement the ``"direct_beam"`` absolute
scaling method.

Method
------

The Master Requirements Document (section 12.2) describes two methods to scale the measured intensity into
absolute units (1/cm): measuring a calibrated standard sample, or measuring the attenuated empty (direct)
beam. The direct beam method is the one most commonly used on GPSANS, and it is the GPSANS default
(``"absoluteScaleMethod": "direct_beam"``).

The empty beam provides a measure of the neutron flux incident on the sample. It is measured with the main
detector in the same instrument configuration as the sample and background, except for two changes:

- the beam stop is moved out of the beam, so that the whole direct beam reaches the detector;
- a calibrated attenuator is placed upstream of the first collimating aperture, so that the pixel count rate
  stays within the linear response range of the detector.

In drtsans the empty beam measurement is the beam center run (``"beamCenter"`` in the reduction parameters).

The summed intensity :math:`I_{BS}` of the beam spot is obtained from the pixels within a radius
:math:`R_{beam}` of the beam center (Eq. 12.5 of the Master Requirements Document):

.. math::

    I_{BS} = \sum_{x,y;\,|\vec{r}_{x,y}| < R_{beam}} I(x,y)

The radius :math:`R_{beam}` is the reduction parameter ``"DBScalingBeamRadius"``, in mm. If it is ``null``,
it is estimated from the source and sample apertures.

The attenuator transmits only a fraction :math:`f(\lambda)` of the neutrons, where :math:`\lambda` is the
wavelength in Å. The intensity of the unattenuated beam :math:`I_T`, and its uncertainty, are:

.. math::

    I_T = \frac{I_{BS}}{f(\lambda)}, \qquad
    \delta I_T = I_T \sqrt{\left(\frac{\delta f}{f}\right)^2 + \left(\frac{\delta I_{BS}}{I_{BS}}\right)^2}

The Master Requirements Document writes this as :math:`I_T = \chi I_{BS}` (Eqs. 12.7 and 12.8), with an
intensity scaling factor :math:`\chi` for the attenuator. The two notations are related by
:math:`\chi = 1/f(\lambda)`.

Finally, the sample intensity, already normalized by the sample thickness, is divided by :math:`I_T`
(Eqs. 12.9 and 12.10), with the relative uncertainties added in quadrature.

Attenuator transmission
-----------------------

The transmitted fraction of each attenuator depends on the wavelength. It is modeled as

.. math::

    f(\lambda) = A e^{-B \lambda} + C

with :math:`\lambda` in Å, :math:`B` in 1/Å, and :math:`A`, :math:`C` dimensionless. The uncertainties
:math:`\delta A`, :math:`\delta B` and :math:`\delta C` of the fitted coefficients are treated as uncorrelated:

.. math::

    \delta f = \sqrt{\left(e^{-B \lambda}\,\delta A\right)^2
    + \left(A \lambda e^{-B \lambda}\,\delta B\right)^2
    + \left(\delta C\right)^2}

The attenuator and the wavelength are read from the ``attenuator`` and ``wavelength`` sample logs of the
empty beam run, not of the sample run. The ``attenuator`` log value identifies the attenuator:

.. list-table::
   :widths: 20 30 50
   :header-rows: 1

   * - Log value
     - Attenuator
     - Transmitted fraction
   * - 0
     - Undefined
     - :math:`f = 1`, :math:`\delta f = 0`
   * - 1
     - Close
     - :math:`f = 1`, :math:`\delta f = 0`
   * - 2
     - Open
     - :math:`f = 1`, :math:`\delta f = 0`
   * - 3
     - x3
     - from the coefficients file
   * - 4
     - x30
     - from the coefficients file
   * - 5
     - x300
     - from the coefficients file
   * - 6
     - x2k
     - from the coefficients file
   * - 7
     - x10k
     - from the coefficients file
   * - 8
     - x100k
     - from the coefficients file
   * - negative
     - Undefined
     - :math:`f = 1`, :math:`\delta f = 0`, with a warning

As a guide to which attenuator a direct beam measurement uses, the Master Requirements Document (section 7.1)
reports that on GPSANS the beam is attenuated by a factor of about 10k or 2k at a wavelength of 4.75 Å. For
wavelengths of 12 Å and longer, it is attenuated by a factor of about 30 with the 40 mm source aperture, and
not attenuated with the 20 mm source aperture.

Attenuation coefficients file
-----------------------------

The coefficients :math:`A`, :math:`B`, :math:`C` and their uncertainties are read from a comma-separated text
file:

- lines that are blank or start with ``#`` are ignored;
- every other line holds the attenuator name (as in the table above) followed by six finite numbers:
  :math:`A`, :math:`\delta A`, :math:`B`, :math:`\delta B`, :math:`C`, :math:`\delta C`;
- each attenuator name appears only once;
- only the attenuators used in the reduction need to be present.

drtsans includes a default file, ``GPSANS_attenuation_coefficients.txt``, with the coefficients provided by
the GPSANS instrument team in 2020:

.. literalinclude:: ../../../src/drtsans/configuration/GPSANS_attenuation_coefficients.txt
   :language: text

The location of the default file in the installed package is given by

.. code-block:: python

    import os
    import drtsans

    print(os.path.join(drtsans.configdir, "GPSANS_attenuation_coefficients.txt"))

A copy of this file is a good starting point for a custom file.

Using a custom coefficients file
--------------------------------

To reduce data with a custom coefficients file, set the reduction parameter
``"AttenuationCoefficientsFileName"`` to the path of the file. A relative path is searched in the current
directory, then in the directories of ``"dataDirectories"``:

.. code-block:: json

    {
      "instrumentName": "GPSANS",
      "configuration": {
        "absoluteScaleMethod": "direct_beam",
        "DBScalingBeamRadius": null,
        "AttenuationCoefficientsFileName": "/path/to/my_attenuation_coefficients.txt"
      }
    }

The attenuation factor of an empty beam workspace can also be computed directly:

.. code-block:: python

    from drtsans.mono.gpsans import attenuation_factor

    factor, factor_error = attenuation_factor(empty_beam_workspace, "/path/to/my_attenuation_coefficients.txt")

The reduction stops with an error if:

- the file does not exist (reported when the reduction parameters are validated);
- a line of the file does not have seven comma-separated fields, its attenuator name is empty or repeated, or one
  of the six coefficients is not a finite number;
- the attenuator of the empty beam run is not listed in the file;
- the ``attenuator`` log value of the empty beam run is positive and not an integer from 0 to 8. This includes
  runs converted from SPICE files with the attenuator open, whose log holds a positive stage position in mm.

With ``"direct_beam"`` scaling, the coefficients file is read before any output file is written, even when the
beam is not attenuated, so an invalid file stops the reduction early.

Attenuation coefficients in the reduction log
---------------------------------------------

When ``"absoluteScaleMethod"`` is ``"direct_beam"``, the reduction log (the ``*_reduction_log.hdf`` file in the
output directory) records the attenuator of the empty beam run and all the coefficients of the attenuation
coefficients file used in the reduction. They are saved next to the absolute scale factor:

.. code-block:: text

    /reduction_information/special_parameters/absolute_scale/
        method                  "direct_beam"
        factor/
            value               absolute scale factor
            error
        attenuation/
            attenuator          attenuator of the empty beam run
            coefficients/
                x3/
                    A/
                        value
                        error
                    B/                  in 1/Å
                        value
                        error
                    C/
                        value
                        error
                x30/
                ...

- ``attenuator`` is the name of the attenuator given by the ``attenuator`` log value of the empty beam run, as in
  the table above. It is ``"Undefined"``, ``"Close"`` or ``"Open"`` when the beam is not attenuated, and
  ``"Undefined"`` for a negative log value.
- ``coefficients`` holds one group for every line of the coefficients file, named after the attenuator,
  including attenuators not used in the reduction. With a custom file, only the attenuators listed in that file
  appear.
- Each coefficient is a group with its ``value`` and ``error``. The HDF5 datasets carry no units: the values
  and errors of :math:`B` are in 1/Å, and :math:`A` and :math:`C` are dimensionless.

No ``attenuation`` group is written when ``"absoluteScaleMethod"`` is ``"standard"``.

For example, to read the attenuator and its coefficients with ``h5py``:

.. code-block:: python

    import h5py

    with h5py.File("/path/to/output/sample_reduction_log.hdf", "r") as log:
        attenuation = log["reduction_information/special_parameters/absolute_scale/attenuation"]
        attenuator = attenuation["attenuator"][()].decode()
        print(f"Attenuator: {attenuator}")
        if attenuator in attenuation["coefficients"]:
            coefficients = attenuation["coefficients"][attenuator]
            for name in ("A", "B", "C"):
                value = coefficients[name]["value"][()]
                error = coefficients[name]["error"][()]
                print(f"{name} = {value} +/- {error}")

Parameters
----------

.. list-table::
   :widths: 25 65 10
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``"absoluteScaleMethod"``
     - Absolute scaling method, ``"standard"`` (multiply by ``"StandardAbsoluteScale"``) or ``"direct_beam"``.
     - ``"direct_beam"``
   * - ``"DBScalingBeamRadius"``
     - Radius :math:`R_{beam}` (mm) of the beam spot summed in the empty beam run. If ``null``, it is
       estimated from the source and sample apertures.
     - ``null``
   * - ``"AttenuationCoefficientsFileName"``
     - Path to the attenuation coefficients file. If ``null``, the default file included in drtsans is used.
     - ``null``
