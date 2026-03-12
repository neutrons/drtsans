"""
Spin state (polarization) event filtering.

This module provides filtering based on polarization device states, splitting events
into scattering cross-sections defined by polarizer and analyzer configurations.
"""

from operator import itemgetter
from typing import List, Optional, Tuple

from mantid.api import AnalysisDataService
from mantid.simpleapi import (
    AddSampleLog,
    CreateEmptyTableWorkspace,
    FilterEvents,
    GenerateEventsFilter,
    GroupWorkspaces,
    logger,
    mtd,
)

from drtsans.filterevents.basefilter import FilterStrategy
from drtsans.polarization import (
    PV_POLARIZER_FLIPPER,
    PV_ANALYZER_FLIPPER,
    PV_POLARIZER,
    PV_POLARIZER_VETO,
    PV_ANALYZER,
    PV_ANALYZER_VETO,
)
from drtsans.samplelogs import SampleLogs
from drtsans.type_hints import MantidWorkspace
from drtsans.dataobjects import workspace_handle


def extract_times(
    times: List[int],
    device_on: bool,
    is_polarizer: Optional[bool] = False,
    is_analyzer: Optional[bool] = False,
    is_polarizer_veto: Optional[bool] = False,
    is_analyzer_veto: Optional[bool] = False,
) -> List[Tuple[int, bool, List[bool]]]:
    """
    Extract time intervals and associate them with specific device states.

    This helper function converts a list of timestamps into a structured format
    that indicates which polarization devices changed state at each time.

    Parameters
    ----------
    times : list of int
        List of timestamps (in nanoseconds) when state changes occurred
    device_on : bool
        Whether these times correspond to the device turning ON (True) or OFF (False)
    is_polarizer : bool, optional
        Whether this state change affects the polarizer. Default is False.
    is_analyzer : bool, optional
        Whether this state change affects the analyzer. Default is False.
    is_polarizer_veto : bool, optional
        Whether this state change affects the polarizer veto. Default is False.
    is_analyzer_veto : bool, optional
        Whether this state change affects the analyzer veto. Default is False.

    Returns
    -------
    list of tuple
        List of tuples, each containing:
        - timestamp (int): Time of the state change in nanoseconds
        - device_on (bool): Whether the device is turning ON or OFF
        - device_mask (list of bool): [is_polarizer, is_analyzer, is_polarizer_veto, is_analyzer_veto]

    Examples
    --------
    >>> times = [1000000, 2000000, 3000000]
    >>> result = extract_times(times, device_on=True, is_polarizer=True)
    >>> # Each timestamp is marked as a polarizer ON event
    >>> assert result[0] == (1000000, True, [True, False, False, False])

    Notes
    -----
    This function is used internally to build a chronological list of all device
    state changes, which is then used to construct the splitter table.
    """
    return [
        (times[i], device_on, [is_polarizer, is_analyzer, is_polarizer_veto, is_analyzer_veto])
        for i in range(len(times))
    ]


def create_table(
    change_list: list, start_time: int, has_polarizer: Optional[bool] = True, has_analyzer: Optional[bool] = True
) -> object:
    """
    Create a Mantid table workspace defining time intervals for cross-sections.

    This function processes a chronological list of device state changes to create
    a splitter table that defines when each cross-section (combination of polarizer
    and analyzer states) is valid.

    Parameters
    ----------
    change_list : list of tuple
        Sorted list of state change tuples from :func:`extract_times`.
        Each tuple contains ``(timestamp, device_on, device_mask)`` — see
        :meth:`SpinFilter._build_change_list` for the full table of possible
        ``device_mask`` combinations and their meanings.
    start_time : int
        Run start time in nanoseconds, used to normalize intervals
    has_polarizer : bool, optional
        Whether a polarizer is present in the experiment. Default is True.
    has_analyzer : bool, optional
        Whether an analyzer is present in the experiment. Default is True.

    Returns
    -------
    Mantid TableWorkspace
        A table with columns:
        - 'start' (float): Start time of the interval in seconds
        - 'stop' (float): Stop time of the interval in seconds
        - 'target' (str): Cross-section label (e.g., 'On_Off', 'Off_On')

    Notes
    -----
    The table only includes time intervals where:
    - Both polarizer and analyzer states are fully specified
    - Neither polarizer nor analyzer veto is active

    Cross-section labels follow the format '{polarizer_state}_{analyzer_state}'
    where each state is either 'On' or 'Off'.

    Time intervals before start_time are discarded or truncated.
    """
    split_table_ws = CreateEmptyTableWorkspace(OutputWorkspace=mtd.unique_hidden_name())
    split_table_ws.addColumn("float", "start")
    split_table_ws.addColumn("float", "stop")
    split_table_ws.addColumn("str", "target")

    # Indices for device_mask
    POLARIZER, ANALYZER, POL_VETO, ANA_VETO = 0, 1, 2, 3

    # Current device states
    current_state = [False] * 4
    current_state_t0 = 0

    # Track whether we have a fully specified state
    # If no device, assume it's specified (don't wait for it)
    specified = [not has_polarizer, not has_analyzer]

    def _add_cross_section_row(start_ns: int, stop_ns: int):
        """Add a row to the table workspace for the current state."""
        xs_label = "%s_%s" % ("On" if current_state[POLARIZER] else "Off", "On" if current_state[ANALYZER] else "Off")
        start = int(start_ns - start_time)
        stop = stop_ns - start_time

        # Discard or truncate intervals before start_time
        if start < 0 and stop <= 0:
            return  # Completely before start
        if start < 0 < stop:
            start = 0.0  # Keep fragment after start_time

        # Mantid expects times in seconds
        split_table_ws.addRow([start * 1e-9, stop * 1e-9, xs_label])

    # Process each state change chronologically
    for change_time, device_on, device_mask in change_list:
        # If we have a complete, valid state, add it to the table
        if (
            specified[POLARIZER]
            and specified[ANALYZER]
            and not current_state[POL_VETO]
            and not current_state[ANA_VETO]
        ):
            _add_cross_section_row(current_state_t0, change_time)

        # Update the current state based on this change
        for i, changed in enumerate(device_mask):
            if changed:
                if i in (POLARIZER, ANALYZER):
                    specified[i] = True
                current_state[i] = device_on

        current_state_t0 = change_time

    return split_table_ws


class SpinFilter(FilterStrategy):
    """
    Filter events by polarization (spin) states.

    This strategy splits event data into scattering cross-sections based on the
    states of polarization devices (polarizer and analyzer). Each cross-section
    corresponds to a specific combination of device states (e.g., both ON, both OFF,
    or one ON and one OFF).

    The filtering is based on sample logs that record the real-time states of the
    polarizer and analyzer flippers, as well as any veto signals that indicate
    periods when the devices are transitioning or unstable.

    Parameters
    ----------
    workspace : str or IEventWorkspace
        The input workspace to be filtered
    pv_polarizer_state : str, optional
        Name of the sample log for polarizer flipper state.
        Default is from drtsans.polarization.PV_POLARIZER_FLIPPER
    pv_analyzer_state : str, optional
        Name of the sample log for analyzer flipper state.
        Default is from drtsans.polarization.PV_ANALYZER_FLIPPER
    pv_polarizer_veto : str, optional
        Name of the sample log for polarizer veto signal.
        Default is from drtsans.polarization.PV_POLARIZER_VETO
    pv_analyzer_veto : str, optional
        Name of the sample log for analyzer veto signal.
        Default is from drtsans.polarization.PV_ANALYZER_VETO
    check_devices : bool, optional
        Whether to check for the presence of polarizer/analyzer in the experiment.
        If False, assumes both devices are present. Default is True.

    Attributes
    ----------
    pv_polarizer_state : str
        Sample log name for polarizer state
    pv_analyzer_state : str
        Sample log name for analyzer state
    pv_polarizer_veto : str
        Sample log name for polarizer veto
    pv_analyzer_veto : str
        Sample log name for analyzer veto
    check_devices : bool
        Whether device presence is checked
    _has_polarizer : bool
        Whether a polarizer was detected in the experiment
    _has_analyzer : bool
        Whether an analyzer was detected in the experiment

    Notes
    -----
    If no polarizer or analyzer is detected, and check_devices is True, the
    raw workspace is returned ungrouped with a warning.

    The filter creates a custom splitter table based on device state transitions
    rather than using Mantid's standard time or log interval filters. This allows
    for complex state combinations and veto handling.

    Cross-section workspaces are labeled as '{polarizer}_{analyzer}' where each
    is either 'On' or 'Off'. For example:
    - 'On_On': Both polarizer and analyzer ON
    - 'On_Off': Polarizer ON, analyzer OFF
    - 'Off_On': Polarizer OFF, analyzer ON
    - 'Off_Off': Both polarizer and analyzer OFF
    """

    def __init__(
        self,
        workspace,
        pv_polarizer_state: str = PV_POLARIZER_FLIPPER,
        pv_analyzer_state: str = PV_ANALYZER_FLIPPER,
        pv_polarizer_veto: str = PV_POLARIZER_VETO,
        pv_analyzer_veto: str = PV_ANALYZER_VETO,
        check_devices: bool = True,
    ):
        """
        Initialize the spin filter.

        Parameters
        ----------
        workspace : str or IEventWorkspace
            The input workspace to filter
        pv_polarizer_state : str, optional
            Sample log name for polarizer flipper state
        pv_analyzer_state : str, optional
            Sample log name for analyzer flipper state
        pv_polarizer_veto : str, optional
            Sample log name for polarizer veto
        pv_analyzer_veto : str, optional
            Sample log name for analyzer veto
        check_devices : bool, optional
            Whether to check for device presence
        """
        super().__init__(workspace)
        self.pv_polarizer_state = pv_polarizer_state
        self.pv_analyzer_state = pv_analyzer_state
        self.pv_polarizer_veto = pv_polarizer_veto
        self.pv_analyzer_veto = pv_analyzer_veto
        self.check_devices = check_devices
        self._has_polarizer = False
        self._has_analyzer = False

    def generate_filter(self) -> Optional[dict]:
        """
        Generate custom splitter table for polarization states.

        This method builds a chronological list of all device state changes and
        creates a custom splitter table defining when each cross-section is valid.

        Returns
        -------
        dict or None
            Empty dict if filtering succeeds (splitter created separately).
            None if no devices are present or no valid intervals exist.

        Notes
        -----
        Unlike other filter strategies, this method creates a custom splitter table
        stored in self.splitter_workspace rather than returning parameters for
        GenerateEventsFilter. The empty dict return value signals apply_filter()
        to skip GenerateEventsFilter and use the custom table directly.
        """
        # Check if devices are present
        sample_logs = SampleLogs(self.workspace)
        if self.check_devices:
            polarizer = sample_logs.get(PV_POLARIZER, 0).value if PV_POLARIZER in sample_logs else 0
            analyzer = sample_logs.get(PV_ANALYZER, 0).value if PV_ANALYZER in sample_logs else 0
        else:
            polarizer = analyzer = 1  # Assume present

        self._has_polarizer = polarizer > 0
        self._has_analyzer = analyzer > 0

        if not self._has_polarizer and not self._has_analyzer:
            logger.warning("No polarizer/analyzer information available")
            return None

        # Build change list from device states
        change_list = self._build_change_list()

        if not change_list:
            return None

        # Create custom splitter table
        start_time = workspace_handle(self.workspace).run().startTime().totalNanoseconds()
        self.splitter_workspace = create_table(
            change_list, start_time, has_polarizer=self._has_polarizer, has_analyzer=self._has_analyzer
        )

        # Return empty dict to signal that custom splitter is ready
        return {}

    def _build_change_list(self) -> List[Tuple[int, bool, List[bool]]]:
        """
        Extract all device state changes from sample logs.

        This method queries the sample logs for polarizer and analyzer state changes,
        as well as veto signals, and combines them into a chronological list.

        Returns
        -------
        list of tuple
            Sorted list of state change tuples, each containing
            ``(timestamp, device_on, device_mask)`` where:

            - ``timestamp`` (int): time of the state change in nanoseconds
            - ``device_on`` (bool): ``True`` if the device turned ON, ``False`` if it turned OFF
            - ``device_mask`` (list of bool): ``[is_polarizer, is_analyzer, is_polarizer_veto, is_analyzer_veto]``.
               Only one of these should be ``True`` for each entry, indicating which device's state changed.

            The possible combinations and their meanings are:

            ============  ============  ============  ==============  ==============  ================================
            ``device_on`` ``is_pol``    ``is_ana``    ``is_pol_veto`` ``is_ana_veto`` Meaning
            ============  ============  ============  ==============  ==============  ================================
            True          True          False         False           False           Polarizer flipper turned **ON**
            False         True          False         False           False           Polarizer flipper turned **OFF**
            True          False         True          False           False           Analyzer flipper turned **ON**
            False         False         True          False           False           Analyzer flipper turned **OFF**
            True          False         False         True            False           Polarizer veto became **active**
            False         False         False         True            False           Polarizer veto **lifted**
            True          False         False         False           True            Analyzer veto became **active**
            False         False         False         False           True            Analyzer veto **lifted**
            ============  ============  ============  ==============  ==============  ================================
        """
        change_list = []

        # Extract polarizer state changes
        if self._has_polarizer:
            change_list.extend(self._extract_device_changes(self.pv_polarizer_state, is_polarizer=True))
            if self.pv_polarizer_veto:
                change_list.extend(self._extract_veto_changes(self.pv_polarizer_veto, is_polarizer_veto=True))

        # Extract analyzer state changes
        if self._has_analyzer:
            change_list.extend(self._extract_device_changes(self.pv_analyzer_state, is_analyzer=True))
            if self.pv_analyzer_veto:
                change_list.extend(self._extract_veto_changes(self.pv_analyzer_veto, is_analyzer_veto=True))

        return sorted(change_list, key=itemgetter(0))

    def _extract_device_changes(self, log_name: str, **kwargs) -> List:
        """
        Extract ON/OFF transitions for a device from its state log.

        Parameters
        ----------
        log_name : str
            Name of the sample log recording device state
        **kwargs
            Device identification flags (is_polarizer or is_analyzer)

        Returns
        -------
        list
            List of state change tuples
        """
        changes = []

        # Device ON (log value ~1.0)
        splitws, _ = GenerateEventsFilter(
            InputWorkspace=self.workspace,
            LogName=log_name,
            MinimumLogValue=0.99,
            MaximumLogValue=1.01,
            TimeTolerance=0,
            UnitOfTime="Nanoseconds",
            OutputWorkspace=mtd.unique_hidden_name(),
            InformationWorkspace=mtd.unique_hidden_name(),
        )
        time_dict = splitws.toDict()
        changes.extend(extract_times(time_dict["start"], True, **kwargs))
        changes.extend(extract_times(time_dict["stop"], False, **kwargs))

        # Device OFF (log value ~0.0)
        splitws, _ = GenerateEventsFilter(
            InputWorkspace=self.workspace,
            LogName=log_name,
            MinimumLogValue=-0.01,
            MaximumLogValue=0.01,
            TimeTolerance=0,
            UnitOfTime="Nanoseconds",
            OutputWorkspace=mtd.unique_hidden_name(),
            InformationWorkspace=mtd.unique_hidden_name(),
        )
        time_dict = splitws.toDict()
        changes.extend(extract_times(time_dict["start"], False, **kwargs))
        changes.extend(extract_times(time_dict["stop"], True, **kwargs))

        return changes

    def _extract_veto_changes(self, log_name: str, **kwargs) -> List:
        """
        Extract veto signal transitions from a veto log.

        Parameters
        ----------
        log_name : str
            Name of the sample log recording veto signal
        **kwargs
            Veto identification flags (is_polarizer_veto or is_analyzer_veto)

        Returns
        -------
        list
            List of veto change tuples
        """
        changes = []

        # Veto ON (log value ~1.0) - vetoed periods to exclude
        splitws, _ = GenerateEventsFilter(
            InputWorkspace=self.workspace,
            LogName=log_name,
            MinimumLogValue=0.99,
            MaximumLogValue=1.01,
            TimeTolerance=0,
            UnitOfTime="Nanoseconds",
            OutputWorkspace=mtd.unique_hidden_name(),
            InformationWorkspace=mtd.unique_hidden_name(),
        )
        time_dict = splitws.toDict()
        changes.extend(extract_times(time_dict["start"], True, **kwargs))
        changes.extend(extract_times(time_dict["stop"], False, **kwargs))

        return changes

    def apply_filter(self, output_workspace: str) -> None:
        """
        Apply polarization filtering with custom splitter table.

        This method overrides the base class implementation to use the custom
        splitter table created in generate_filter() rather than using
        GenerateEventsFilter.

        Parameters
        ----------
        output_workspace : str
            Name for the output workspace group
        """
        filter_params = self.generate_filter()

        # If no filtering possible, return raw workspace
        if filter_params is None:
            GroupWorkspaces([self.workspace], OutputWorkspace=output_workspace)
            return

        # Check if splitter table has content
        if self.splitter_workspace.rowCount() == 0:
            logger.warning("No valid cross-section intervals found")
            GroupWorkspaces([self.workspace], OutputWorkspace=output_workspace)
            return

        # Use custom splitter table
        correction_workspace = mtd.unique_hidden_name()
        outputs = FilterEvents(
            InputWorkspace=self.workspace,
            SplitterWorkspace=self.splitter_workspace,
            GroupWorkspaces=True,
            FilterByPulseTime=False,
            OutputWorkspaceIndexedFrom1=False,
            CorrectionToSample="None",
            SpectrumWithoutDetector="Skip",
            SplitSampleLogs=True,
            RelativeTime=True,
            ExcludeSpecifiedLogs=True,
            OutputTOFCorrectionWorkspace=correction_workspace,
            OutputWorkspaceBaseName=output_workspace,
        )
        AnalysisDataService.remove(correction_workspace)

        # Add cross-section IDs to each workspace
        for ws in outputs[-1]:
            xs_id = str(ws).replace(f"{output_workspace}_", "")
            AddSampleLog(Workspace=ws, LogName="cross_section_id", LogText=xs_id)
            # These two lines mimic what FilterEvents does when passed an information workspace
            ws.setComment(xs_id)
            ws.setTitle(xs_id)

    def inject_metadata(self, workspace: MantidWorkspace) -> None:
        """
        Inject metadata into all polarization-filtered cross-sections.

        Adds common metadata (slice number, total slices) and polarization-specific
        metadata to each cross-section workspace in the group.

        Parameters
        ----------
        workspace : MantidWorkspace
            The workspace group (or its name) containing the filtered cross-sections
        """
        workspace_group = workspace_handle(workspace)
        num_slices = workspace_group.getNumberOfEntries()

        for n in range(num_slices):
            slice_workspace = workspace_group.getItem(n)
            samplelogs = SampleLogs(slice_workspace)

            # Add common metadata
            samplelogs.insert("slice", n + 1)
            samplelogs.insert("number_of_slices", num_slices)

            # Add polarization-specific metadata
            samplelogs.insert("slice_parameter", "polarization_state")

            # Extract cross-section ID from workspace
            xs_id = (
                samplelogs.get("cross_section_id", "Unknown").value if "cross_section_id" in samplelogs else "Unknown"
            )
            samplelogs.insert("cross_section", xs_id)

            # Add device presence information
            samplelogs.insert("has_polarizer", int(self._has_polarizer))
            samplelogs.insert("has_analyzer", int(self._has_analyzer))
