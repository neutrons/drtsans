.. _user.corrections.transmission:

Transmission
============

The measured intensity must be corrected to account for sample transmission:

.. math::

    I'_{sample}(x,y,\lambda) = \frac{I_{sample}(x,y,\lambda)}{T(x,y,\lambda)}

The transmission correction is calculated using a transmission run and an empty beam transmission
run (reference).
It is possible to specify the transmission value directly in the transmission parameter ``"value"``,
however, this is mostly used for diagnostic purposes.

.. code-block:: json

    {
      "sample": {
        "runNumber": 10010,
        "loadOptions": {},
        "thickness": 1.0,
        "transmission": {
          "runNumber": null,
          "value": null,
          "errorTolerance": null
        }
      },
      "background": {
        "runNumber": null,
        "transmission": {
          "runNumber": null,
          "value": null
        }
      },
      "emptyTransmission": {
        "runNumber": null,
        "value": null
      },
      "mmRadiusForTransmission": null,
      "useTimeSliceTransmission": false,
      "useThetaDepTransCorrection": true
    }

Time Slice Transmission (Bio-SANS only)
---------------------------------------

When using time slicing (``"useTimeSlice": true``), users can optionally calculate the transmission
correction using the time sliced sample run by setting the parameter
``"useTimeSliceTransmission": true``. The sample transmission run number is ignored when
``"useTimeSliceTransmission": true``. The time slice transmission option can be used when the sample
transmission is expected to change over time.

Time slices for which the transmission calculation fails will be skipped. The transmission
calculation can fail due to all transmission values being NaN or if the transmission error is
higher than the allowed relative transmission error, which is configurable in the sample
transmission parameter ``"errorTolerance"`` (default: 0.01). For example, the last time slice may
be shorter and, therefore, include fewer counts, resulting in large statistical errors in the
transmission calculation.

.. code-block:: json

    {
      "sample": {
        "runNumber": 10010,
        "loadOptions": {},
        "thickness": 1.0,
        "transmission": {
          "runNumber": null,
          "value": null,
          "errorTolerance": 0.05
        }
      },
      "emptyTransmission": {
        "runNumber": 10005,
        "value": null
      },
      "useTimeSlice": true,
      "timeSliceInterval": 100.0,
      "useTimeSliceTransmission": true
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
     - Beam radius within which the transmission will be calculated. If ``null``, then the beam
       radius is calculated from the sample logs.
     - ``null``
   * - ``"useThetaDepTransCorrection"``
     - If ``true``, a theta dependent transmission correction will be applied, which takes into
       account the effect of the scattering angle on the transmission.
     - ``true``
   * - ``"useTimeSliceTransmission"``
     - (`Only for Bio-SANS and when` ``"useTimeSlice": true``.) If ``true``, the transmission
       correction will be calculated using the time sliced sample run itself instead of a separate
       sample transmission run. This is useful when the sample transmission is expected to change
       over time. Slices with relative transmission error larger than
       ``"transmissionErrorTolerance"`` will be skipped.
     - ``false``
   * - ``"errorTolerance"``
     - (`Only for Bio-SANS and when` ``"useTimeSlice": true`` `and` ``"useTimeSliceTransmission": true``.)
       Maximum relative transmission error.
     - 0.01
