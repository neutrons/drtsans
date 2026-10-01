.. _user.corrections.attenuation:

Absolute Scaling with an Attenuated Direct Beam (GPSANS)
========================================================

This correction applies only to GPSANS. BIOSANS and EQSANS do not implement the ``"direct_beam"`` absolute
scaling method.

.. |br| raw:: html

   <br />

.. topic:: On this page

   - `Methods`_ |br|
     How the attenuated direct beam is used for absolute scaling.
   - `Attenuator Transmission`_ |br|
     How the attenuation factor and its uncertainty are computed.
   - `Attenuation Coefficients File`_ |br|
     Format of the default calibration-history file and custom override files.
   - `Using a Custom Coefficients File`_ |br|
     How to point a GPSANS reduction at user-supplied attenuation coefficients.
   - `Attenuation Coefficients in the Reduction Log`_ |br|
     Where the selected formula and coefficients are recorded in reduction output.
   - `Parameters`_ |br|
     Reduction parameters that control attenuated direct-beam scaling.

Methods
-------

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

Attenuator Transmission
-----------------------

The transmitted fraction of each attenuator may depend on the wavelength. The model is read from the attenuation
coefficients file. The current default calibration uses

.. math::

    f(\lambda) = A e^{-B \lambda} + C

with :math:`\lambda` in Å, :math:`B` in 1/Å, and :math:`A`, :math:`C` dimensionless. More generally, the file
contains a formula using fitted parameter names and the optional reserved variable ``wavelength``. If
``wavelength`` is omitted, the attenuation is wavelength independent.

The uncertainties of the fitted parameters are treated as uncorrelated and propagated automatically:

.. math::

    \delta f = \sqrt{\sum_p \left(\frac{\partial f}{\partial p}\delta p\right)^2}

where :math:`p` runs over the fitted parameters in the formula. No wavelength uncertainty is included.

The attenuator and, when needed, the wavelength are read from the ``attenuator`` and ``wavelength`` sample logs
of the empty beam run, not of the sample run. The ``attenuator`` log value identifies the attenuator:

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
   * - between 0 and 3, not an integer
     - Undefined
     - :math:`f = 1`, :math:`\delta f = 0`, with a warning

As a guide to which attenuator a direct beam measurement uses, the Master Requirements Document (section 7.1)
reports that on GPSANS the beam is attenuated by a factor of about 10k or 2k at a wavelength of 4.75 Å. For
wavelengths of 12 Å and longer, it is attenuated by a factor of about 30 with the 40 mm source aperture, and
not attenuated with the 20 mm source aperture.

Attenuation Coefficients File
-----------------------------

The formula and fitted parameter values are read from a text file:

- lines that are blank or start with ``#`` are ignored;
- ``formula =`` declares the SymPy-style attenuation formula;
- parameter lines use ``parameter = value, uncertainty``;
- the first attenuator block declares the fitted parameter names, and every later attenuator block must define
  the same parameter names;
- supported formula syntax is numbers, fitted parameter names, ``wavelength``, arithmetic operators,
  parentheses, and the functions ``abs``, ``acos``, ``asin``, ``atan``, ``cos``, ``cosh``, ``erf``, ``exp``,
  ``log``, ``log10``, ``sin``, ``sinh``, ``sqrt``, ``tan`` and ``tanh``;
- the constant ``pi`` is supported;
- each attenuator name appears only once.

drtsans includes a timestamped default file, ``GPSANS_attenuation_coefficients.txt``. Each calibration block
starts with ``[effective YYYY-MM-DD]``. During reduction, the block is selected from the empty beam timestamp,
using ``start_time`` first, then ``run_start``, then ``run_begin``. The most recent block whose effective date is
not later than the run date is used. The existing packaged coefficients are effective ``1990-01-01``:


.. literalinclude:: ../../../src/drtsans/configuration/GPSANS_attenuation_coefficients.txt
   :language: text

The location of the default file in the installed package is given by

.. code-block:: python

    import os
    import drtsans

    print(os.path.join(drtsans.configdir, "GPSANS_attenuation_coefficients.txt"))

A copy of the parameter blocks in this file is a good starting point for a custom file. Custom files do not use
``[effective ...]`` sections, because they are explicit overrides for the reduction.

Using a Custom Coefficients File
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
- the custom file uses the old seven-column comma-separated format;
- the file has malformed formula, attenuator, or parameter lines;
- the attenuator blocks do not all define the same fitted parameter names;
- a value or uncertainty is not a finite number;
- all default-file run timestamp logs, ``start_time``, ``run_start`` and ``run_begin``, are missing, or the
  selected timestamp is earlier than the earliest effective calibration block;
- the attenuator of the empty beam run is not listed in the file;
- the ``attenuator`` log value of the empty beam run is 3 or more and not an integer from 3 to 8. This includes
  runs converted from SPICE files with the attenuator open, whose log holds a positive stage position in mm.

With ``"direct_beam"`` scaling, the coefficients file is read before any reduced I(Q) output is written, even when
the beam is not attenuated, so an invalid file stops the reduction early.

Attenuation Coefficients in the Reduction Log
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
            fit_function        attenuation formula
            effective_date      selected default block date, omitted for custom files
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
  ``"Undefined"`` for a negative log value or a non-integer log value between 0 and 3.
- ``coefficients`` holds one group for every line of the coefficients file, named after the attenuator,
  including attenuators not used in the reduction. With a custom file, only the attenuators listed in that file
  appear.
- Each coefficient is a group with its ``value`` and ``error``. Parameter names come from the selected formula
  block. The HDF5 datasets carry no units; interpret units from the formula.

No ``attenuation`` group is written when ``"absoluteScaleMethod"`` is ``"standard"``.

For example, to read the attenuator and its coefficients with ``h5py``:

.. code-block:: python

    import h5py

    with h5py.File("/path/to/output/sample_reduction_log.hdf", "r") as log:
        attenuation = log["reduction_information/special_parameters/absolute_scale/attenuation"]
        attenuator = attenuation["attenuator"][()].decode()
        formula = attenuation["fit_function"][()].decode()
        print(f"Attenuator: {attenuator}")
        print(f"Formula: {formula}")
        if attenuator in attenuation["coefficients"]:
            coefficients = attenuation["coefficients"][attenuator]
            for name in coefficients:
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
