"""
Base strategy class and factory for event filtering.

This module defines the abstract base class for all event filtering strategies,
along with utility functions for determining which strategy to use based on
reduction configuration parameters.
"""

from abc import ABC, abstractmethod
from typing import Generator, Optional, Union, Tuple
from mantid.api import IEventWorkspace
from mantid.simpleapi import GenerateEventsFilter, FilterEvents

from drtsans.dataobjects import workspace_handle
from drtsans.samplelogs import SampleLogs
from drtsans.type_hints import MantidWorkspace


class FilterStrategy(ABC):
    """
    Abstract base class for event filtering strategies.

    This class defines the interface that all concrete filtering strategies must implement.
    It provides a common framework for generating filters, applying them to workspaces,
    and extracting metadata from the filtered results.

    Parameters
    ----------
    workspace : str or IEventWorkspace
        The input workspace to be filtered

    Attributes
    ----------
    workspace : str or IEventWorkspace
        The input workspace being filtered
    splitter_workspace : str
        Name of the temporary workspace holding the splitter table
    info_workspace : str
        Name of the temporary workspace holding filter information
    """

    FILTER_WORKSPACE_NAME = "_filter"
    INFO_WORKSPACE_NAME = "_info"

    def __init__(self, workspace: Union[str, IEventWorkspace]):
        """
        Initialize the filter strategy.

        Parameters
        ----------
        workspace : str or IEventWorkspace
            The input workspace to be filtered
        """
        self.workspace: str = str(workspace)  # store the name
        self.splitter_workspace: str = self.FILTER_WORKSPACE_NAME
        self.info_workspace: str = self.INFO_WORKSPACE_NAME

    @abstractmethod
    def generate_filter(self) -> Optional[dict]:
        """
        Generate filter parameters for Mantid's GenerateEventsFilter algorithm.

        This method must be implemented by concrete strategy classes to provide
        the specific parameters needed for their filtering approach.

        Returns
        -------
        dict or None
            Parameters to pass to GenerateEventsFilter algorithm.
            Return None if filtering is not applicable (e.g., no polarization devices found).
        """
        pass

    @abstractmethod
    def inject_metadata(self, workspace: MantidWorkspace) -> None:
        """
        Inject metadata into all sliced workspaces in the output workspace group.

        This method must be implemented by concrete strategy classes to add
        filter-specific metadata to each workspace slice. It should add both
        common metadata (slice number, total slices) and filter-specific metadata
        (parameter name, interval, start/end values with appropriate units).

        Parameters
        ----------
        workspace : MantidWorkspace
            The workspace group (or its name) containing the filtered workspaces

        Notes
        -----
        Each implementation should add the following metadata to each slice:
        - 'slice': Slice number (1-based)
        - 'number_of_slices': Total number of slices
        - 'slice_info': Comment from the workspace
        - 'slice_parameter': Name/description of the slicing parameter
        - 'slice_interval': The interval value used
        - 'slice_start': Start value for this slice (with units if applicable)
        - 'slice_end': End value for this slice (with units if applicable)
        """
        pass

    def _inject_common_metadata(
        self, workspace: MantidWorkspace
    ) -> Generator[Tuple[int, SampleLogs, str], None, None]:
        """
        Inject common metadata into each slice and yield per-slice objects for further use.

        Handles the shared metadata block — slice number, total slice count, and workspace
        comment — for every workspace in the group, then yields ``(n, samplelogs, slice_info)``
        so the calling ``inject_metadata`` override can add filter-specific entries without
        repeating the common block.

        Parameters
        ----------
        workspace : MantidWorkspace
            The workspace group (or its name) containing the filtered workspaces

        Yields
        ------
        n : int
            Zero-based index of the current slice within the group.
        samplelogs : SampleLogs
            Sample-log accessor for the current slice workspace, already populated with
            ``"slice"``, ``"number_of_slices"``, and ``"slice_info"``.
        slice_info : str
            Raw comment string from the current slice workspace, as returned by
            ``getComment()``.
        """
        workspace_group = workspace_handle(workspace)
        num_slices = workspace_group.getNumberOfEntries()
        for n in range(num_slices):
            slice_workspace = workspace_group.getItem(n)
            samplelogs = SampleLogs(slice_workspace)
            samplelogs.insert("slice", n + 1)
            samplelogs.insert("number_of_slices", num_slices)
            slice_info = slice_workspace.getComment()
            samplelogs.insert("slice_info", slice_info)
            yield n, samplelogs, slice_info

    def apply_filter(self, output_workspace: str) -> None:
        """
        Apply the filtering logic to create split workspaces.

        This method can be overridden by subclasses that need custom filtering
        behavior (e.g., polarization filtering with custom splitter tables).
        The default implementation uses Mantid's standard GenerateEventsFilter
        and FilterEvents algorithms.

        Parameters
        ----------
        output_workspace : str
            Name for the output workspace group that will contain filtered workspaces

        Notes
        -----
        After filtering, empty workspaces should be removed from the group.
        Subclasses overriding this method should maintain similar behavior.
        """
        filter_params = self.generate_filter()

        # Allow subclasses to handle no-filter case
        if filter_params is None:
            return

        GenerateEventsFilter(
            InputWorkspace=self.workspace,
            OutputWorkspace=self.splitter_workspace,
            InformationWorkspace=self.info_workspace,
            **filter_params,
        )

        FilterEvents(
            InputWorkspace=self.workspace,
            OutputWorkspaceBaseName=output_workspace,
            SplitterWorkspace=self.splitter_workspace,
            InformationWorkspace=self.info_workspace,
            FilterByPulseTime=True,
            GroupWorkspaces=True,
            OutputWorkspaceIndexedFrom1=True,
        )
