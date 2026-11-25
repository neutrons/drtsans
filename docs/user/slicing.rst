.. _user.slicing:

Slicing
=======

The user can perform slicing of a data set by either time or sample log (but not both). Data
reduction is performed separately on each slice according to the specified reduction parameters.
The output files will have a slice index appended to the configured ``"outputFileName"``, e.g.
``myFileName_0_Iq.dat``, ``myFileName_1_Iq.dat``, etc.

For EQSANS run in frame-skip mode, each frame will have its own time slices.
Thus, output files will also contain the frame index, e.g.
``myFileName_0_frame_0_Iq.dat`` and ``myFileName_0_frame_1_Iq.dat``.
If wedge binning is selected, the wedge label appears between the slice and frame labels,
e.g. ``myFileName_0_wedge_0_frame_0_Iq.dat``. PNG plot filenames follow the
same convention as the corresponding reduced data files.

.. contents:: Contents
   :local:

Time Slicing
------------

The user can specify time slicing to slice a run by time intervals and, optionally, a time period.
The following parameters specify the time slicing:

.. code-block:: json

    {
        "useTimeSlice": true,
        "timeSliceInterval": 4800,
        "timeSliceOffset": 0.0,
        "timeSlicePeriod": null,
    }

``timeSliceInterval`` is either a single slice length or a list of slice lengths, in seconds (see Example 3).
``timeSliceOffset`` is the start of the first slice, in seconds after the run start.

The final slice
+++++++++++++++

Without ``timeSlicePeriod``, the slices are consecutive and cover the run from ``timeSliceOffset`` to the run end.
A single number ``t`` is the same as the list ``[t]``. The values of the list are used in order, and the last value
``t`` is repeated until the run is covered. The time left at the run end, shorter than ``t``, is not made into a
slice of its own: the final slice absorbs it, so the final slice lasts at least ``t`` and less than ``2t``.

When the run is shorter than the values of the list, the run end cuts the current slice short, however short that
is, and the remaining values are not used. An error message in the log reports it.

The final slice ends one pulse period after the last pulse recorded in the data file, when that is later than the
run end, so that it keeps the events of the last pulse. The pulses are those of the data acquisition, 1/60 s apart:
at SNS (EQSANS) they are the neutron pulses of the source, and at HFIR (BIOSANS and GPSANS), a reactor, they are
written by the data acquisition. The final slice is therefore a fraction of a second longer than the rule above
gives: for example 90.065 s instead of 90.048 s for 70 s slices of a 300.048 s run.


Example 1: Slicing using time interval
++++++++++++++++++++++++++++++++++++++

.. code-block:: json

    {
        "useTimeSlice": true,
        "timeSliceInterval": 2.0,
        "timeSliceOffset": 0.0,
        "timeSlicePeriod": null,
    }

.. figure:: media/time_slicing_interval.png
   :alt: Diagram of time slicing using an interval
   :width: 800px

Example 2: Slicing using time interval and offset
+++++++++++++++++++++++++++++++++++++++++++++++++

.. code-block:: json

    {
        "useTimeSlice": true,
        "timeSliceInterval": 2.0,
        "timeSliceOffset": 1.0,
        "timeSlicePeriod": null,
    }

.. figure:: media/time_slicing_offset.png
   :alt: Diagram of time slicing using an interval and offset
   :width: 800px

Events recorded before ``timeSliceOffset`` (here 0-1 s) are discarded. The remaining 9 s are not an exact multiple
of ``timeSliceInterval``: instead of a 1 s slice at the run end, the final slice absorbs the leftover and lasts 3 s,
from 7 s to 10 s (see `The final slice`_).

Example 3: Slicing using a list of slice lengths
++++++++++++++++++++++++++++++++++++++++++++++++

``timeSliceInterval`` can be a list of slice lengths, in seconds. Each value gives one slice, in order, starting at
``timeSliceOffset``.

.. code-block:: json

    {
        "useTimeSlice": true,
        "timeSliceInterval": [2.0, 4.0],
        "timeSliceOffset": 0.0,
        "timeSlicePeriod": null,
    }

.. figure:: media/time_slicing_interval_list.png
   :alt: Diagram of time slicing using a list of slice lengths
   :width: 800px

The 13 s run is longer than the 6 s of the list, so the last value, 4 s, is repeated, and the final slice absorbs
the 3 s left at the run end: it lasts 7 s, from 6 s to 13 s. A warning in the log reports the slices added to the
list. For a 5 s run, the slices would be 0-2 s and 2-5 s, the run end cutting the second slice short.

The offset delays only the first slice and leaves the lengths unchanged: with ``"timeSliceOffset": 1.0``, the 13 s
run would give the slices 1-3 s, 3-7 s and 7-13 s, the final slice absorbing 2 s.

A list of more than one value cannot be combined with ``timeSlicePeriod``: the reduction parameters are rejected
when they are validated.

Example 4: Slicing using time slice period
++++++++++++++++++++++++++++++++++++++++++

When using a time slice period, slices for a given time interval are reduced together.
In the figure, red slices are reduced together, and blue slices are reduced together.

.. code-block:: json

    {
        "useTimeSlice": true,
        "timeSliceInterval": 2.0,
        "timeSliceOffset": 0.0,
        "timeSlicePeriod": 4.0,
    }

.. figure:: media/time_slicing_period.png
   :alt: Diagram of time slicing using an interval and a period
   :width: 800px

The offset also applies to periodic slicing. With ``"timeSliceOffset": 1.0`` and the same interval and
period, the first period starts at 1 s and the last window (9-10 s) is shorter:

.. code-block:: json

    {
        "useTimeSlice": true,
        "timeSliceInterval": 2.0,
        "timeSliceOffset": 1.0,
        "timeSlicePeriod": 4.0,
    }

.. figure:: media/time_slicing_offset_period.png
   :alt: Diagram of time slicing using an interval, an offset and a period
   :width: 800px

Note: if the period is not an integer multiple of the time slice interval, then the remainder of each period
forms its last, shorter slice: `The final slice`_ rule does not apply within a period. In the example below, each
period holds two 3 second slices, then one 2 second slice.

.. code-block:: json

    {
        "useTimeSlice": true,
        "timeSliceInterval": 3.0,
        "timeSliceOffset": 0.0,
        "timeSlicePeriod": 8.0,
    }

Example 5: Time slicing and event filtering during loading
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The time slicing parameters are applied relative to the whole data from the beginning to the end of
the run, independent of any event filtering during loading (i.e. using the parameter
``"loadOptions"``).

This example shows EQSANS (SNS time-of-flight instrument) configuration:

.. code-block:: json

    {
        "sample": {
            "loadOptions": {"FilterByTimeStart": 2.0, "FilterByTimeStop": 8.71}
        },
        "configuration": {
            "useTimeSlice": true,
            "timeSliceInterval": 2.0,
            "timeSliceOffset": 0.0,
            "timeSlicePeriod": 4.0,
        }
    }

.. figure:: media/time_slicing_period_and_event_filtering.png
   :alt: Diagram of log slicing
   :width: 800px

Event Filtering for BIOSANS and GPSANS
++++++++++++++++++++++++++++++++++++++++

BIOSANS and GPSANS (HFIR monochromatic instruments) also support event filtering using
``loadOptions``. The same ``FilterByTimeStart`` and ``FilterByTimeStop`` parameters can be used
to filter events before reduction.

Example for BIOSANS or GPSANS:

.. code-block:: json

    {
        "instrumentName": "BIOSANS",
        "sample": {
            "runNumber": "1322",
            "loadOptions": {"FilterByTimeStart": 10.0, "FilterByTimeStop": 300.0}
        }
    }

As for EQSANS (Example 5), the load-time cutoffs and the time slices are independent. Events outside
``FilterByTimeStart`` and ``FilterByTimeStop`` are discarded during loading, but the slices still start at
``timeSliceOffset`` after the run start, whatever the value of ``FilterByTimeStart``. Slices that fall
outside the cutoffs hold no events; see `Empty and partially filled slices`_.

.. note::
   The ``FilterByTimeStart`` and ``FilterByTimeStop`` parameters specify times in seconds relative
   to the start of the run. When the run is split (time slicing, log slicing, or polarization on GPSANS), these
   options are passed to ``LoadEventNexus`` for all instruments. When the run is not split, they are passed
   to ``LoadEventAsWorkspace2D`` for HFIR instruments (BIOSANS and GPSANS) and to ``LoadEventNexus`` for
   EQSANS.

Empty and partially filled slices
+++++++++++++++++++++++++++++++++

Slices that hold no events, because they fall outside the load-time cutoffs or during a period without
beam, are dropped and not reduced. The output suffix (``_0``, ``_1``, ...) counts only the remaining
slices, so ``_0`` is not necessarily the first time window of the run. If only one slice remains, the
output files have no suffix.

A slice can also be skipped during the reduction, with a warning in the log: on BIOSANS and GPSANS when its
monitor counts are zero, on BIOSANS also when its time-slice transmission fails the transmission checks, and
on EQSANS when its neutron flux is zero. The suffixes are not renumbered, so the written files then have a
gap in their suffixes where a slice was skipped.

A slice that straddles a cutoff keeps only part of its window. Its ``duration`` and ``proton_charge``
sample logs cover only the time it actually holds, so the normalization uses its actual exposure.

In the example below, a 60 s run is sliced every 10 s and only the events after 25 s are loaded. The windows
0-10 s and 10-20 s are dropped. The window 20-30 s, which holds the events from 25 s to 30 s, is written
with suffix ``_0``, and the windows 30-40 s, 40-50 s and 50-60 s with suffixes ``_1``, ``_2`` and ``_3``.

.. code-block:: json

    {
        "sample": {
            "loadOptions": {"FilterByTimeStart": 25.0}
        },
        "configuration": {
            "useTimeSlice": true,
            "timeSliceInterval": 10.0,
            "timeSliceOffset": 0.0,
            "timeSlicePeriod": null,
        }
    }

.. figure:: media/time_slicing_empty_slices.png
   :alt: Diagram of time slicing with empty slices removed
   :width: 800px

With periodic slicing, an empty window does not remove a slice, because each slice also collects the events
of the same window in the other periods. In the figure of Example 5, the empty window 0-2 s belongs to slice 0,
which still holds the events recorded from 4 s to 6 s and from 8 s to 8.71 s.

Slice metadata
++++++++++++++

Each slice receives the following sample logs:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Sample log
     - Meaning
   * - ``slice_start``, ``slice_end``
     - Time slicing: start and end of the slice's window, in seconds from the run start. This is the full
       window even when the slice is only partly filled, and the final slice may end slightly after the run end
       (see `The final slice`_).
       Periodic slicing: start and end of the window, in seconds from the start of each period. Log slicing: lower and upper
       values of the sample log, in the log's units.
   * - ``slice_interval``
     - Time slicing: the length of the slice, in seconds; it differs between slices for a list of slice lengths
       and for the final slice, whose length includes the fraction of a second after the run end (see
       `The final slice`_). Periodic slicing: the value of ``timeSliceInterval``, in seconds. Log slicing: the
       log-value interval.
   * - ``slice_period``
     - The value of ``timeSlicePeriod``, in seconds (periodic slicing only).
   * - ``slice_parameter``
     - ``relative time from start`` (time slicing), ``periodic time from start`` (periodic slicing), or the
       name of the sample log (log slicing).
   * - ``slice``, ``number_of_slices``
     - Position of the slice, counting from 1, and number of slices. Both count only the remaining slices.
   * - ``duration``
     - The exposure the slice actually covers, in seconds.

To place a slice in time, use ``slice_start`` and ``slice_end``, not the output suffix.

Log Slicing
-----------

Log slicing refers to slicing by the value of a `sample log
<https://docs.mantidproject.org/nightly/tutorials/mantid_basic_course/connecting_data_to_instruments/06_sample_logs.html>`_
in the NeXus file corresponding to the sample run.

.. code-block:: json

    {
        "useLogSlice": true,
        "logSliceName": "MySampleLog",
        "logSliceInterval": 200.0,
    }

Example: Log slicing using log "Temperature"
++++++++++++++++++++++++++++++++++++++++++++

.. code-block:: json

    {
        "useLogSlice": true,
        "logSliceName": "Temperature",
        "logSliceInterval": 10.0,
    }

.. figure:: media/log_slicing_interval.png
   :alt: Diagram of log slicing
   :width: 800px


Slicing and Error-Weighted Binning
----------------------------------

Note that using error-weighted averaging in combination with slicing may result in disparity
between time sliced and non-time sliced data if there is not enough data in the slices for the
error distribution to converge.

.. code-block:: json

    {
        "useErrorWeighting": true,
    }

Figure 1 compares I(Q) for the whole run (blue dot markers) to I(Q) for time slices
of the data (solid lines) for the time slice intervals 600 s, 1200 s, 2400 s and 4800 s. The time
sliced I(Q) converges to the total I(Q) with increasing time slice interval.

.. figure:: media/AlCell_He_1bar_timeslice_vs_no_timeslice.png
   :alt: I(Q) with and without time slicing for different time slice intervals using error-weighted averaging.
   :width: 800px

   Figure 1: I(Q) for time slices of a run compared to the total run for different time slice intervals and
   with error-weighted averaging (``"useErrorWeighting": true``).

In Figure 2, the same comparison between time sliced and non-time sliced data but with
error-weighted averaging turned off (i.e. using arithmetic averaging) shows no dependence on the
time slice interval.

.. figure:: media/AlCell_He_1bar_timeslice_vs_no_timeslice_no_error_weighting.png
   :alt: I(Q) with and without time slicing for different time slice intervals using arithmetic averaging.
   :width: 800px

   Figure 2: I(Q) for time slices of a run compared to the total run for different time slice intervals and
   with error-weighted averaging (``"useErrorWeighting": false``).
