.. _user.corrections.transmission:

Transmission
============

The measured intensity must be corrected to account for sample transmission:

.. math::

    I'_{sample}(x,y,\lambda) = \frac{I_{sample}(x,y,\lambda)}{T(\lambda,x,y)}

The transmission correction is calculated using a transmission run and an empty beam transmission
run (reference).
It is possible to specify the transmission value directly in the transmission parameter ``"value"``,
however, this is mostly used for diagnostic purposes.

.. code-block:: json

    {
      "sample": {
        "runNumber": None,
        "loadOptions": {},
        "thickness": 1.0,
        "transmission": {
          "runNumber": None,
          "value": None,
          "errorTolerance": None
        }
      },
      "background": {
        "runNumber": None,
        "transmission": {
          "runNumber": None,
          "value": None
        }
      },
      "emptyTransmission": {
        "runNumber": None,
        "value": None
      },
      "mmRadiusForTransmission": None,
      "useTimeSliceTransmission": False,
      "useThetaDepTransCorrection": True,
    }

Time Slice Transmission (Bio-SANS only)
---------------------------------------

When using time slicing (``"useTimeSlice": True``), users can optionally calculate the transmission
correction using the time sliced sample run by setting the parameter
``"useTimeSliceTransmission": True``. The sample transmission run number is ignored when
``"useTimeSliceTransmission": True``. The time slice transmission option can be used when the sample
transmission is expected to change over time.

Time slices for which the transmission calculation fails will be skipped. The transmission
calculation can fail due to all transmission values being NaN or if the transmission error is
higher than the allowed relative transmission error, which is configurable in the sample
transmission parameter ``"errorTolerance"`` (default: 0.01).

.. code-block:: json

    {
      "sample": {
        "runNumber": 10010,
        "loadOptions": {},
        "thickness": 1.0,
        "transmission": {
          "runNumber": None,
          "value": None,
          "errorTolerance": 0.05
        }
      },
      "emptyTransmission": {
        "runNumber": 10005,
        "value": None
      },
      "useTimeSlice": True,
      "timeSliceInterval": 100.0,
      "useTimeSliceTransmission": True,
    }

Parameters
----------

.. list-table::
   :widths: 25 65 10
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``"mmRadiusForTransmission"``
     - Beam radius within which the transmission will be calculated. If ``None``, then the beam
       radius is calculated from the sample logs.
     - ``None``
   * - ``"useThetaDepTransCorrection"``
     - If ``True``, a theta dependent transmission correction will be applied, which takes into
       account the effect of the scattering angle on the transmission.
     - ``True``
   * - ``"useTimeSliceTransmission"``
     - (`Only for Bio-SANS and when` ``"useTimeSlice": True``.) If ``True``, the transmission
       correction will be calculated using the time sliced sample run itself instead of a separate
       sample transmission run. This is useful when the sample transmission is expected to change
       over time. Slices with relative transmission error larger than
       ``"transmissionErrorTolerance"`` will be skipped.
     - ``False``
   * - ``"errorTolerance"``
     - (`Only for Bio-SANS and when` ``"useTimeSlice": True`` `and` ``"useTimeSliceTransmission": True``.) Maximum relative transmission error.
     - 0.01
