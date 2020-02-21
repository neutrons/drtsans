import copy
import enum
import itertools
import json
import numpy as np
import numexpr
import os
import sys

r""" Hyperlinks to mantid algorithms
ApplyCalibration <https://docs.mantidproject.org/nightly/algorithms/ApplyCalibration-v1.html>
CloneWorkspace <https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html>
CreateEmptyTableWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateEmptyTableWorkspace-v1.html>
CreateWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>
DeleteWorkspaces <https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspaces-v1.html>
FilterEvents <https://docs.mantidproject.org/nightly/algorithms/FilterEvents-v1.html>
GenerateEventsFilter <https://docs.mantidproject.org/nightly/algorithms/GenerateEventsFilter-v1.html>
Integration <https://docs.mantidproject.org/nightly/algorithms/Integration-v1.html>
Load <https://docs.mantidproject.org/nightly/algorithms/Load-v1.html>
LoadEmptyInstrument <https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html>
LoadInstrument <https://docs.mantidproject.org/nightly/algorithms/LoadInstrument-v1.html>
LoadNexus <https://docs.mantidproject.org/nightly/algorithms/LoadNexus-v1.html>
MaskDetectors <https://docs.mantidproject.org/nightly/algorithms/MaskDetectors-v1.html>
MaskDetectorsIf <https://docs.mantidproject.org/nightly/algorithms/MaskDetectorsIf-v1.html>
ReplaceSpecialValues <https://docs.mantidproject.org/nightly/algorithms/ReplaceSpecialValues-v1.html>
SaveNexus <https://docs.mantidproject.org/nightly/algorithms/SaveNexus-v1.html>
"""
from mantid.simpleapi import (ApplyCalibration, CloneWorkspace, CreateEmptyTableWorkspace, CreateWorkspace,
                              DeleteWorkspaces, FilterEvents, GenerateEventsFilter, Integration, Load,
                              LoadEmptyInstrument, LoadInstrument, LoadNexus, MaskDetectors, MaskDetectorsIf,
                              ReplaceSpecialValues, SaveNexus)
from mantid.api import mtd

r"""
Hyperlinks to drtsans functions
namedtuplefy, unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
SampleLogs <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
TubeCollection <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tubecollection.py>
"""  # noqa: E501
from drtsans.instruments import InstrumentEnumName, instrument_enum_name, instrument_standard_name
from drtsans.path import exists as file_exists
from drtsans.settings import namedtuplefy, unique_workspace_dundername
from drtsans.samplelogs import SampleLogs
from drtsans.tubecollection import TubeCollection


__all__ = ['apply_calibrations', ]

r"""Flags a problem when running the barscan algorithm that identifies the pixel corresponding
to the bottom of the shadow cast by the bar on the detector array."""
INCORRECT_PIXEL_ASSIGNMENT = -1

r"""Default files storing the metadata of the pixel calibrations. There's one file for each instrument."""
database_file = {InstrumentEnumName.BIOSANS: '/HFIR/CG3/shared/calibration/pixel_calibration.json',
                 InstrumentEnumName.EQSANS: '/SNS/EQSANS/shared/calibration/pixel_calibration.json',
                 InstrumentEnumName.GPSANS: '/HFIR/CG2/shared/calibration/pixel_calibration.json'}


class CalType(enum.Enum):
    r"""Enumerate the possible types of pixel calibrations"""
    BARSCAN = 'BARSCAN'
    TUBEWIDTH = 'TUBEWIDTH'


def day_stamp(input_workspace):
    r"""
    Find the day stamp (e.g 20200311 for March 11, 2020) using the "start_time" metadata from the
    Nexus events file as input.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        Workspace from which the day stamp is to be retrieved.

    Returns
    -------
    int
    """
    return int(SampleLogs(input_workspace).start_time.value[0:10].replace('-', ''))


class CalibrationNotFound(Exception):
    """Exception to be raised when no appropriate calibration is found in the database"""
    pass


class Table:
    r"""Container for a table of pixel calibration item data, plus metadata

    The table object holds two attributes:
    - metadata, dict, informs about the calibration run, instrument, detector array.
    - table, ~mantid.api.TableWorkspace, containing the actual calibration data.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    metadata: dict
        Dictionary with the following fields about the calibration:
        - caltype, str, the type of calibration (BARSCAN, TUBEWIDTH)
        - instrument, str, standard name of the instrument for which the calibration was carried out.
        - component, str, standard name of the double detector array for which the calibration was carried out.
        - daystamp, int, 8-digit integer whose digits are to be understood as YYYYMMDD.
        - run_numbers, list, list of run numbers that encompassed the calibration.
    detector_ids: list
        List of detector IDs for which a calibration has been carried out.
    positions: list
        List of Y-coordinates for each detector, in meters.
    heights: list
        List of detector heights (along the Y-), in meters.
    widths: list
        List of detector widths (along the X-axis), in meters.
    """

    @classmethod
    def compose_table_name(cls, metadata):
        r"""Standard workspace name for a calibration table, built as a composite name using the
        calibration type, instrument, component, and daystamp. (e.g. "barscan_gpsans_detector1_20200311")

        Parameters
        ----------
        metadata: dict
            Dictionary containing the metadata of one calibration

        Returns
        -------
        str
        """
        m = metadata  # handy shortcut
        return f'{m["caltype"].lower()}_{m["instrument"]}_{m["component"]}_{str(m["daystamp"])}'

    @classmethod
    def load(cls, database, caltype, instrument, component, daystamp, output_workspace=None):
        r"""
        Load a Nexus file containing a calibration into a ```Table``` object.

        **Mantid algorithms used:**
        :ref:`LoadNexus <algm-LoadNexus-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/LoadNexus-v1.html>

        Parameters
        ----------
        database: str
            Path to JSON file containing metadata for different calibrations.
        caltype: str
            Type of calibration (BARSCAN, TUBEWIDHT).
        instrument: str
            Standard name of the instrument for which the calibration was carried out.
        component: str
            Standard name of the double detector array for which the calibration was carried out.
        daystamp: int
            8-digit integer whose digits are to be understood as YYYYMMDD. The returned calibration
            will have a daystamp equal or more recent.
        output_workspace: str
            Name of the output ~mantid.api.TableWorkspace containing the calibration data. If
            :py:obj:`None`, a composite name is created using the calibration type, instrument, component,
            and daystamp. (e.g. "barscan_gpsans_detector1_20200311").

        Returns
        -------
        ~drtsans.pixel_calibration.Table
        """
        # Search the database for a match to the required metadata
        not_found_message = f'No suitable {caltype}_{instrument}_{component} calibration found in {database}'
        with open(database, mode='r') as json_file:
            entries = json.load(json_file)  # list of metadata entries
            required = {caltype, instrument, component}  # plausible metadata entries must match these metadata pieces
            candidates = [entry for entry in entries if required.issubset(set(entry.keys()))]
            if len(candidates) == 0:
                raise CalibrationNotFound(not_found_message)
            candidates.sort(key=lambda c: c['daystamp'])  # sort candidates by increasing daystamp
        for i, candidate in enumerate(candidates):
            if candidate['daystamp'] > daystamp:
                if i == 0:
                    raise CalibrationNotFound(not_found_message)
                metadata = candidates[i - 1]
                if output_workspace is None:
                    output_workspace = Table.compose_table_name(metadata)
                table = LoadNexus(metadata['tablefile'], OutputWorkspace=output_workspace)
                return Table(table, metadata)
        raise CalibrationNotFound(not_found_message)

    @classmethod
    def build_mantid_table(cls, output_workspace, detector_ids, positions=None, heights=None, widths=None):
        r"""
        Instantiate a Table workspace with input calibration data.

        **Mantid algorithms used:**
        :ref:`CreateEmptyTableWorkspace <algm-CreateEmptyTableWorkspace-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/CreateEmptyTableWorkspace-v1.html>

        Parameters
        ----------
        output_workspace: str
            Name of the output table workspace.
        detector_ids: list
            List of detector IDs for which a calibration has been carried out.
        positions: list
            List of Y-coordinates for each detector, in meters.
        heights: list
            List of detector heights (along the Y-), in meters.
        widths: list
            List of detector widths (along the X-axis), in meters.

        Returns
        -------
        ~mantid.api.TableWorkspace
        """
        columns_data = {'Detector Y Coordinate': positions, 'Detector Height': heights, 'Detector Width': widths}
        [columns_data.pop(column) for column in list(columns_data.keys()) if columns_data[column] is None]
        table = CreateEmptyTableWorkspace(OutputWorkspace=output_workspace)
        table.addColumn(type='int', name='Detector ID')
        [table.addColumn(type='double', name=column) for column in columns_data]
        for i in range(len(detector_ids)):
            row = {'Detector ID': detector_ids[i]}
            row.update({column: data[i] for column, data in columns_data.items()})
            table.addRow(row)
        return table

    @classmethod
    def validate_metadata(cls, metadata):
        r"""
        Verify the metadata contains entries for the instrument, double-detector-array, and day stamp.

        Parameters
        ----------
        metadata: dict

        Returns
        -------
        bool

        Raises
        ------
        ValueError
            The metadata is missing one of the required entries.
        """
        required_keys = {'instrument', 'component', 'daystamp'}
        if required_keys.issubset(metadata.keys()) is False:
            raise ValueError(f'Metadata is missing one or more of these entries: {required_keys}')

    def __init__(self, metadata, detector_ids, positions=None, heights=None, widths=None):
        Table.validate_metadata(metadata)
        output_workspace = Table.compose_table_name(metadata)
        self.table = Table.build_mantid_table(output_workspace, detector_ids,
                                              positions=positions, heights=heights, widths=widths)
        self.metadata = copy.copy(metadata)

    def __getattr__(self, item):
        r"""Serve metadata's keys as attributes of the ```Table``` object"""
        if item not in self.__dict__:
            return self.__dict__['metadata'][item]
        return self.__dict__[item]

    def column_values(self, name):
        r"""
        Return a list of values for the selected table column.

        Possible names are 'Detector ID', 'Detector Y Coordinate', 'Detector Height', and 'Detector Width'.

        Parameters
        ----------
        name: str
            Name of the column. Must match the name of one of the columns in the ~mantid.api.TableWorkspace
            ```table``` attribute.

        Returns
        -------
        list
        """
        column_names = self.table.getColumnNames()
        try:
            column_index = column_names.index(name)
        except ValueError:
            raise ValueError(f'"{name}" is not a column name of the calibration table')
        return self.table.column(column_index)

    @property
    def detector_ids(self):
        r"""List of pixel positions stored in the calibration table."""
        return self.column_values('Detector ID')

    @property
    def positions(self):
        r"""List of pixel positions stored in the calibration table."""
        return self.column_values('Detector Y Coordinate')

    @property
    def heights(self):
        r"""List of pixel heights stored in the calibration table."""
        return self.column_values('Detector Height')

    @property
    def widths(self):
        r"""List of pixel widths stored in the calibration table."""
        return self.column_values('Detector Width')

    def apply(self, input_workspace, output_workspace=None):
        r"""
        Apply a calibration to an input workspace and return the calibrated workspace.

        **Mantid algorithms used:**
        :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html>
        :ref:`ApplyCalibration <algm-ApplyCalibration-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/ApplyCalibration-v1.html>

        Parameters
        ----------
        input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
            Workspace to which calibration needs to be applied.
        output_workspace: str
            Name of the output workspace with calibrated pixels. If :py:obj:`None`, the pixels
            of the input workspace will be calibrated.

        Returns
        -------
        ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        """
        if output_workspace is None:
            output_workspace = str(input_workspace)
        else:
            CloneWorkspace(InputWorkspace=input_workspace, OutputWorkspace=output_workspace)
        ApplyCalibration(Workspace=output_workspace, CalibrationTable=self.table)
        return mtd[output_workspace]

    def save(self, database=None, tablefile=None):
        r"""
        Save the metadata in a JSON file and the table workspace in a Nexus file.

        **Mantid algorithms used:**
        :ref:`SaveNexus <algm-SaveNexus-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/SaveNexus-v1.html>

        Parameters
        ----------
        database: str
            Path to the JSON file where the ```metadata``` dictionary will be appended. If :py:obj:`None`,
            then the appropriate default file from ~drtsans.pixel_calibration.database_file is used.
        tablefile: str
            Path to the Nexus file storing the pixel calibration data. If :py:obj:`None`, then
            a composite name is created using the calibration type, instrument, component,
            and daystamp. (e.g. "barscan_gpsans_detector1_20200311"). The file is saved under
            subdirectory 'calibrations', located within the directory of the ```database``` file.
        """
        if database is None:
            database = database_file[instrument_enum_name(self.instrument)]
        if tablefile is None:
            cal_dir = os.path.join(os.path.dirname(database), 'calibrations')
            if os.path.isdir(cal_dir) is False:
                os.mkdir(cal_dir)
            tablefile = os.path.join(cal_dir, Table.compose_table_name(self.metadata))
        self.metadata['tablefile'] = tablefile
        SaveNexus(InputWorkspace=self.table, Filename=tablefile)
        entries = list()
        if os.path.exists(database):
            with open(database, mode='r') as json_file:
                entries = json.load(json_file)  # list of metadata entries
        entries.append(self.metadata)
        with open(database, mode='w') as json_file:
            json.dump(entries, json_file)

    def as_intensities(self):
        r"""
        Creates one workspace for each pixel property that is calibrated, and the calibration datum is
        saved as the value of the intensity for that pixel. Useful to visualize the calibration in
        MantidPlot's instrument viewer.

        For example, a BARSCAN calibration will generate workspaces ```tablename_positions```
        and ```tablename_heights```, where ```tablename``` is the name of the ~mantid.api.TableWorkspace
        holding the calibration data.

        **Mantid algorithms used:**
        :ref:`CreateWorkspace <algm-CreateWorkspace-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>
        :ref:`LoadEmptyInstrument <algm-LoadEmptyInstrument-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html>
        :ref:`LoadInstrument <algm-LoadInstrument-v1>`,
        <https://docs.mantidproject.org/nightly/algorithms/LoadInstrument-v1.html>
        """
        empty_instrument = LoadEmptyInstrument(InstrumentName=self.instrument,
                                               OutputWorkspace=unique_workspace_dundername())
        detector_ids = self.detector_ids
        calibration_properties = ['positions', 'heights'] if self.caltype == 'BARSCAN' else ['widths', ]
        for cal_prop in calibration_properties:
            values = getattr(self, cal_prop)
            intensities = empty_instrument.extractY()  # extractY returns a copy
            for workspace_index in range(empty_instrument.getNumberHistograms()):
                detector = empty_instrument.getDetector(workspace_index)
                try:
                    row_index = detector_ids.index(detector.getID())
                except ValueError:  # This detector was not calibrated, thus is not in detector_ids
                    continue
                intensities[workspace_index] = values[row_index]
            workspace = CreateWorkspace(DataX=[0, 1], DataY=intensities, Nspec=empty_instrument.getNumberHistograms(),
                                        OutputWorkspace=f'{self.table.name()}_{cal_prop}')
            LoadInstrument(Workspace=workspace, InstrumentName=self.instrument, RewriteSpectraMap=True)
        empty_instrument.delete()


def load_calibration(input_workspace, caltype, component='detector1', database=None, output_workspace=None):
    r"""
    Load a calibration into a ~drtsans.pixel_calibration.Table object.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        Workspace from which calibration session is to be retrieved.
    caltype: str
        Either 'BARSCAN' or 'TUBEWIDTH'.
    component: str
        Name of one of the double detector array panels.
    database: str
        Path to database file containing the metadata for the calibrations. If :py:obj:`None`, the default database
        is used.
    output_workspace: str
        Name of the table workspace containing the calibration session values. If :py:obj:`None`, then a composite
        name is created using the calibration type, instrument, component, and daystamp. (e.g.
        "barscan_gpsans_detector1_20200311")

    Returns
    -------
    ~drtsans.pixel_calibration.Table
    """
    enum_instrument = instrument_enum_name(input_workspace)
    if database is None:
        database = database_file[enum_instrument]
    return Table.load(database, caltype, str(enum_instrument), component, day_stamp(input_workspace),
                      output_workspace=output_workspace)


def apply_calibrations(input_workspace, database=None, calibrations=[cal.name for cal in CalType],
                       output_workspace=None):
    r"""
    Load and apply one or more calibrations to an input workspace.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
        Input workspace whose pixels are to be calibrated.
    database: str
        Path to JSON file containing metadata for different past calibrations.
    calibrations: str, list
        One or more of 'BARSCAN' and/or 'TUBEWIDTH'.
    output_workspace: str
         Name of the output workspace with calibrated pixels. If :py:obj:`None`, the pixels
        of the input workspace will be calibrated.

    Returns
    -------
        ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    if isinstance(calibrations, str):  # we passed only one calibration
        calibrations = [calibrations, ]
    components = {InstrumentEnumName.BIOSANS: ['detector1', 'wing_detector'],
                  InstrumentEnumName.EQSANS: ['detector1'],
                  InstrumentEnumName.GPSANS: ['detector1']}
    for caltype in calibrations:
        for component in components:
            try:
                calibration = load_calibration(input_workspace, caltype, component, database=database)
                calibration.apply(input_workspace)
            except CalibrationNotFound as e:
                sys.stderr.write(e)
    return mtd[output_workspace]


def _consecutive_true_values(values, how_many, reverse=False, raise_message=None):
    r"""
    Find first array index of consecutive `how_many` True values.

    devs - Andrei Savici <saviciat@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    values: list
        list of `True` and `False` items
    how_many: int
        Number of desired consecutive `True` values
    raise_message: str
        Exception message. No exception if :py:obj:`None`, but INCORRECT_PIXEL_ASSIGNMENT is returned

    Returns
    -------
    int

    Raises
    ------
    IndexError
        If no index is found for any of the edges
    RuntimeError
        If a faulty tube is found
    """
    # use the array or the reverse one
    truth_array = values[::-1] if reverse else values
    # create a sub-array of length how_many of True values that we want to find
    pattern = [True]*how_many
    # loop over the input data and return the first index where the next
    # how_many elements match the pattern
    for i in range(len(truth_array) - how_many):
        if truth_array[i: i + how_many] == pattern:
            return len(values) - i - 1 if reverse else i
    # raise an error if the pattern is not found
    else:
        if raise_message is not None:
            raise IndexError(raise_message)
        return INCORRECT_PIXEL_ASSIGNMENT  # signal for non-identified value


@namedtuplefy
def find_edges(intensities, tube_threshold=0.2, shadow_threshold=0.3,
               tube_edge_min_width=3, shadow_edge_min_width=4,
               min_illuminated_length=7):
    r"""
    Find the active length of the tube and the shadow region

    All pixel indexes start from the bottom of the tube, with the first
    index being zero.

    devs - Andrei Savici <saviciat@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    intensities: list
        pixel pixel_intensities along the tube.
    tube_threshold: float
        fraction of the average intensity to determine the tube edges.
    shadow_threshold: float
        fraction of the average intensity to determine the shadow region.
    tube_edge_min_width: int
        required minimum number of consecutive good pixels above
        the tube threshold
    shadow_edge_min_width: int
        required minimum number of consecutive shadowed pixels
    min_illuminated_length: int
        minimum number of illuminated pixels on the active length

    Returns
    -------
    namedtuple
        the fields of the name tuple are:
        - bottom_pixel: first illuminated pixel
        - top_pixel: last illuminated pixel
        - bottom_shadow_pixel: first shadowed pixel
        - above_shadow_pixel= first illuminated pixel above the shadow region
    """
    # calculate minimum intensity thresholds for tube ends and shadows
    average_intensity = np.average(intensities)
    end_threshold = tube_threshold * average_intensity
    shadow_threshold = shadow_threshold * average_intensity

    # Find edges of the tube: want at least tube_edge_min_width pixels
    # (starting from the top or bottom of a tube) that have pixel_intensities greater
    # than the threshold

    illuminated = [bool(i > end_threshold) for i in intensities]
    # The bottom pixel is the first illuminated pixel. It is required that the next tube_edge_min_width pixels are
    # also illuminated
    bottom_pixel = _consecutive_true_values(illuminated, tube_edge_min_width,
                                            raise_message='Could not find bottom tube edge')
    # The top pixel is the last illuminated pixel. It is required that the previous tube_edge_min_width pixels are
    # also illuminated.
    top_pixel = _consecutive_true_values(illuminated, tube_edge_min_width,
                                         raise_message='Could not find top tube edge', reverse=True)

    # Find the shadow region: similar to tube edges, but in this case
    # we want shadow_edge_min_width pixel_intensities less than the shadow threshold,
    # followed by at least one intensity greater than the threshold

    # The bottom pixel shadowed by the bar is the first pixel below the intensity threshold. We require that the
    # next shadow_edge_min_width are also shadowed.
    shadowed = [bool(i < shadow_threshold) for i in intensities[bottom_pixel:]]
    bottom_shadow_pixel = bottom_pixel +\
        _consecutive_true_values(shadowed, shadow_edge_min_width, raise_message='Could not find bottom shadow edge')

    # Find the first illuminated pixel above the bar.
    illuminated = [bool(i > shadow_threshold) for i in
                   intensities[bottom_shadow_pixel + shadow_edge_min_width: top_pixel + 1]]
    # Don't raise if the pixel is not found
    above_shadow_pixel = bottom_shadow_pixel + shadow_edge_min_width +\
        _consecutive_true_values(illuminated, 1, raise_message=None)

    # Check for a faulty tube: we want a certain number of pixels not in the bar shaddow
    active_tube_length = top_pixel - bottom_pixel + 1
    shadow_length = above_shadow_pixel - bottom_shadow_pixel
    if active_tube_length < min_illuminated_length + shadow_length:
        raise RuntimeError('Faulty tube found')

    return dict(bottom_pixel=bottom_pixel, top_pixel=top_pixel,
                bottom_shadow_pixel=bottom_shadow_pixel,
                above_shadow_pixel=above_shadow_pixel)


@namedtuplefy
def fit_positions(edge_pixels, bar_positions, tube_pixels=256, order=5, ignore_value=INCORRECT_PIXEL_ASSIGNMENT):
    r"""
    Fit the position and heights of the pixels in a tube. The bar_positions as a function of
    edge pixels are fitted to a nth order polynomial (by default n=5). The positions of the pixels along the
    tube are the values of the polynomial at integer points, while the heights are the derivatives.

    Description from the master requirements document, section A2.1

    All pixel indexes start from the bottom of the tube, with the first
    index being zero.

    Uses :ref:`~numpy.polynomial.polynomial.polyfit`.

    devs - Andrei Savici <saviciat@ornl.gov>,

    Parameters
    ----------
    edge_pixels: list (or numpy array)
        the bottom pixel for each bar position, as found in `find_edges` function
    bar_positions: list (or numpy array)
        the bar position from the logs for each file in the bar scan
    tube_pixels: integer
        number of pixels for which to calculate positions and heights
    order: integer
        the order of polynomial to be used in the fit (default 5)
    ignore_value: int
        certain positions of the bar (close to the top and bottom of the tube) results in incorrect assignment of the
        edge pixel. In those cases it is expected that the edge pixel has a particular value that flags incorrect
        assignment. The default value is INCORRECT_PIXEL_ASSIGNMENT. These edge pixels will be
        ignored when carrying out the fit.

    Returns
    -------
    namedtuple
        the fields of the name tuple are:
        - calculated_positions: calculated positions of the pixels
        - calculated_heights: calculated pixel heights
    """
    message_len = 'The positions of the bar and edge pixels have to be the same length'
    assert len(edge_pixels) == len(bar_positions), message_len

    # Ignore the incorrectly assigned edge pixels
    edge_pixels = np.array(edge_pixels)
    bar_positions = np.array(bar_positions)
    valid_edge_pixels = edge_pixels[np.where(edge_pixels != ignore_value)]
    valid_bar_positions = bar_positions[np.where(edge_pixels != ignore_value)]

    try:
        # fit the bar positions to a 5th degree polynomial in edge_pixels
        coefficients = np.polynomial.polynomial.polyfit(valid_edge_pixels, valid_bar_positions, int(order))
        # calculate the coefficients of the derivative
        deriv_coefficients = np.polynomial.polynomial.polyder(coefficients)
        # evalutae the positions
        calculated_positions = np.polynomial.polynomial.polyval(np.arange(tube_pixels), coefficients)
        # evaluate the heights
        calculated_heights = np.polynomial.polynomial.polyval(np.arange(tube_pixels), deriv_coefficients)
    except Exception:
        calculated_positions = np.ones(tube_pixels) * np.nan
        calculated_heights = np.ones(tube_pixels) * np.nan

    return dict(calculated_positions=calculated_positions, calculated_heights=calculated_heights,
                coefficients=coefficients)


def event_splitter(barscan_file, split_workspace=None, info_workspace=None, bar_position_log='dcal_Readback'):
    r"""
    Split a Nexus events file containing a full bar scan.

    It is assumed that the bar is shifted by a fixed amount every time we go on to the next scan.

    Parameters
    ----------
    barscan_file: str
        Path to barscan run file containing multiple positions of the bar.
    split_workspace: str
        Name of the table workspace to be used as event splitter workpsce in algorithm FilterEvents. If
        :py:obj:`None`, a random name will be provided.
    info_workspace: str
        Name of the table workspace to be used along ``split_workspace`` in algorithm FilterEvents. If
        :py:obj:`None`, a random name will be provided.
    bar_position_log: str
        Name of the log entry in the barscan run file containing the position of the bar (Y-coordinate, in 'mm')
        with respect to some particular frame of reference, not necessarily the one located at the sample.
    Returns
    -------
    list
        List of bar positions
    """
    if split_workspace is None:
        split_workspace = unique_workspace_dundername()
    if info_workspace is None:
        info_workspace = unique_workspace_dundername()

    barscans_workspace = unique_workspace_dundername()
    Load(barscan_file, OutputWorkspace=barscans_workspace)

    # Find the amount by which the position of bar is shifted every time we go on to the next scan.
    bar_positions = SampleLogs(barscans_workspace)[bar_position_log].value
    bar_delta_positions = bar_positions[1:] - bar_positions[:-1]  # list of shifts in the position of the bar
    bar_delta_positions = bar_delta_positions[bar_delta_positions > 0]  # only shifts where the bar position increases
    # the most likely shift of the bar position, within two significant figures.
    bar_step = float(np.bincount(np.round(100 * bar_delta_positions).astype('int')).argmax()) / 100.

    GenerateEventsFilter(InputWorkspace=barscans_workspace, OutputWorkspace=split_workspace,
                         InformationWorkspace=info_workspace, UnitOfTime='Nanoseconds',
                         LogName=bar_position_log, LogValueInterval=bar_step, LogValueTolerance=bar_step / 2,
                         MinimumLogValue=min(bar_positions), MaximumLogValue=max(bar_positions))

    # Read bar positions from info_workspace using the second column. That column has entries as strings of the form:
    # "Log.dcal_Readback.From.{min}.To.{max}.Value-change-direction:both"
    bar_positions = list()
    for min_max_bar_position in mtd[info_workspace].column(1):
        min_bar_position = float(min_max_bar_position.split('.From.')[-1].split('.To.')[0])
        max_bar_position = float(min_max_bar_position.split('.To.')[-1].split('.Value')[0])
        bar_positions.append((min_bar_position + max_bar_position) / 2)

    return bar_positions


def barscan_workspace_generator(barscan_files, bar_position_log='dcal_Readback'):
    r"""
    A generator to be used in an iteration over the runs that held the bar at a fixed position, returning at each
    iteration the name of the workspace containing a single run (the set of runs make up the bar scan).

    Parameters
    ----------
     barscan_files: str, list
        Path(s) to barscan run file(s). If only one file, it should contain multiple positions of the bar.
        If a list of files, then each file contains the pixel_intensities recorded with a constant position for the
        bar.
    bar_position_log: str
        Name of the log entry in the barscan run file containing the position of the bar (Y-coordinate, in 'mm')
        with respect to some particular frame of reference, not necessarily the one located at the sample.

    Returns
    -------
    tuple
        A two-item tuple containing, in this order:
        - the position of the bar as stated in the logs the run currently being returned.
        - the name of the workspace containing the run currently being returned.
    """
    temporary_workspaces = list()  # store the names of the workspaces to be removed

    def temporary_workspace():
        r"""Random workspace name, and flag it for removal"""
        name = unique_workspace_dundername()
        temporary_workspaces.append(name)
        return name

    if isinstance(barscan_files, str):  # the whole barscan is contained in a single file. Must be splitted
        spliter_workspace = temporary_workspace()
        info_workspace = temporary_workspace()
        # Table to split the barscan run into subruns, each with a fixed position of the bar
        bar_positions = event_splitter(barscan_files, split_workspace=spliter_workspace,
                                       info_workspace=info_workspace, bar_position_log=bar_position_log)
        barscans_workspace = temporary_workspace()
        Load(barscan_files, OutputWorkspace=barscans_workspace)
        splitted_workspace_group = unique_workspace_dundername()
        FilterEvents(InputWorkspace=barscans_workspace,
                     SplitterWorkspace=spliter_workspace, InformationWorkspace=info_workspace,
                     OutputWorkspaceBaseName=splitted_workspace_group,  GroupWorkspaces=True,
                     TimeSeriesPropertyLogs=[bar_position_log], ExcludeSpecifiedLogs=False)
        temporary_workspaces.append(splitted_workspace_group)
        temporary_workspaces.append('TOFCorrectWS')  # spurious workspace spawned by FilterEvents
        for i, bar_position in enumerate(bar_positions):
            yield bar_position, splitted_workspace_group + '_' + str(i)
    else:
        barscan_workspace = temporary_workspace()
        for file_name in barscan_files:
            Load(file_name, OutputWorkspace=barscan_workspace)
            bar_position = SampleLogs(barscan_workspace).find_log_with_units(bar_position_log, 'mm')
            yield bar_position, barscan_workspace
    DeleteWorkspaces(temporary_workspaces)  # clean up the now useless workspaces


def calculate_barscan_calibration(barscan_files, component='detector1', bar_position_log='dcal_Readback',
                                  formula='565 - {y}', order=5):
    r"""
    Calculate pixel positions (only Y-coordinae) as well as pixel heights from a barscan calibration session.

    **Mantid Algorithms used:**
    :ref:`Load <algm-Load-v1>`,

    devs - Andrei Savici <saviciat@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    barscan_files: str, list
        Path(s) to barscan run file(s). If only one file, it should contain multiple positions of the bar.
        If a list of files, then each file contains the pixel_intensities for every pixel for a fixed position of the
        bar.
    component: str
        Name of the detector panel scanned with the bar. Usually, 'detector1`.
    bar_position_log: str
        Name of the log entry in the barscan run file containing the position of the bar (Y-coordinate, in 'mm')
        with respect to some particular frame of reference, not necessarily the one located at the sample.
    formula: str
        Formula to obtain the position of the bar (Y-coordinate) in the frame of reference located at the sample.
    order: int
        Highest degree for the polynomial that will fit the observed positions of the bar.

    Returns
    -------
    dict
        Dictionary with the following entries:
        - instrument, str, Standard name of the instrument.
        - component, str, name of the double detector array, usually "detector1".
        - run, int, run number associated to the calibration.
        - unit: str, the units for the positions and heights. Set to 'mm' for mili-meters.
        - positions, list, List of Y-coordinate for each pixel.
        - heights, list, List of pixel heights.
    """
    instrument_name, number_pixels_in_tube, number_tubes = None, None, None
    run_numbers, daystamp = [], None
    bar_positions = []  # Y-coordinates of the bar for each scan
    # 2D array defining the position of the bar on the detector, in pixel coordinates
    # The first index corresponds to the Y-axis (along each tube), the second to the X-axis (across tubes)
    # Thus, bottom_shadow_pixels[:, 0] indicates bottom shadow pixel coordinates along the very first tubes
    bottom_shadow_pixels = []
    for bar_position, barscan_workspace in barscan_workspace_generator(barscan_files,
                                                                       bar_position_log=bar_position_log):
        if instrument_name is None:
            instrument_name = instrument_standard_name(barscan_workspace)
            daystamp = day_stamp(barscan_workspace)
        run_numbers.append(int(SampleLogs(barscan_workspace).single_value('run_number')))
        # Find out the Y-coordinates of the bar in the reference-of-frame located at the sample
        formula_bar_position_inserted = formula.format(y=bar_position)
        bar_positions.append(float(numexpr.evaluate(formula_bar_position_inserted)))
        bottom_shadow_pixels_per_scan = []  # For the current scan, we have one bottom shadow pixel for each tube
        # A TubeCollection is a list of TubeSpectrum objects, representing a physical tube. Here we obtain the
        # list of tubes for the main double-detector-panel.
        # The view 'decreasing X' sort the tubes by decreasing value of their corresponding X-coordinate. In this view,
        # a double detector panel looks like a single detector panel. When looking at the panel standing at the
        # sample, the leftmost tube has the highest X-coordinate, so the 'decreasing X' view orders the tubes
        # from left to right.
        if number_pixels_in_tube is None:
            # We create a tube collection to figure out the pixel indexes for each tube
            collection = TubeCollection(barscan_workspace, component).sorted(view='decreasing X')
            # pixel_indexes is a list of length equal the number of tubes. Each list element is a list containing
            # the pixel indexes for a particular tube.
            pixel_indexes = [tube.spectrum_info_index for tube in collection]
            number_tubes, number_pixels_in_tube = len(collection), len(collection[0])
            detector_ids = list(itertools.chain.from_iterable(tube.detector_ids for tube in collection))
        pixel_intensities = np.sum(mtd[barscan_workspace].extractY(), axis=1)  # integrated intensity on each pixel
        for pixel_indexes_in_tube in pixel_indexes:  # iterate over each tube, retrieving its pixel indexes
            try:
                # Find the bottom shadow pixel for the current tube and current barscan run
                pixel_intensities_in_tube = pixel_intensities[pixel_indexes_in_tube]
                bottom_shadow_pixels_per_scan.append(find_edges(pixel_intensities_in_tube).bottom_shadow_pixel)
            except Exception:
                # this tube may be malfunctioning for the current barscan
                bottom_shadow_pixels_per_scan.append(INCORRECT_PIXEL_ASSIGNMENT)
        bottom_shadow_pixels.append(bottom_shadow_pixels_per_scan)
    bottom_shadow_pixels = np.array(bottom_shadow_pixels)
    bar_positions = np.array(bar_positions)

    # Deal with corner cases not resolved with the find_edges algorithm
    resolve_incorrect_pixel_assignments(bottom_shadow_pixels, bar_positions)

    if len(bottom_shadow_pixels) <= order:
        raise ValueError(f"There are not enough bar positions to fo a fit with a polynomyal of order {order}.")

    # fit pixel positions for each tube
    positions, heights = [], []
    for tube_index in range(number_tubes):  # iterate over the tubes in the collection
        # Fit the pixel numbers and Y-coordinates of the bar for the current tube with a polynomial
        fit_results = fit_positions(bottom_shadow_pixels[:, tube_index], bar_positions, order=order,
                                    tube_pixels=number_pixels_in_tube)
        # Store the fitted Y-coordinates and heights of each pixel in the current tube
        # Store as lists so that they can be easily serializable
        positions.extend(list(1.e-03 * fit_results.calculated_positions))  # store with units of meters
        heights.extend(list(1.e-03 * fit_results.calculated_heights))  # store with units of meters

    metadata = dict(caltype=CalType.BARSCAN.name,
                    instrument=instrument_name,
                    component=component,
                    daystamp=daystamp,
                    runnumbers=sorted(run_numbers))
    return Table(metadata, detector_ids, positions=positions, heights=heights)


def resolve_incorrect_pixel_assignments(bottom_shadow_pixels, bar_positions):
    r"""
    # Corner case 1:
    # When the bar approaches the bottom of the tubes, the bottom edge of its shadow blends with the non-functioning
    # pixels at the bottom of the tubes. In this scenario, the bottom edge of the bar is assigned as the first
    # non-functioning pixel at the top of the tube. Thus, we must find out for each tube a sudden jump in the
    # identified bottom shadow pixel.
    # Example: the bottom shadow pixels for the first tube as we change the position of the bar have been identified:
    #    249, 230, 209, 187, 165, 145, 129, 103, 82, 67, 52, 39, 21, 249, 249
    # In the example, the last two bottom pixels are erroneous. They must be set to -1
    #
    # Corner case 2:
    # When the identified bottom pixel in a tube is far from the boundary between the illuminated pixels and the non
    # illuminated pixels, it leaves a few non-illuminated pixels as active pixels. The bottom shadow pixel can then
    # be incorrectly assigned as one of these non-illuminated pixels. The results are outliers pixel indexes.
    # Example: 249 248 246 4 244 242.... Here "4" is in incorrectly assigned pixel index as bottom of the bar
    # We must find these outliers and assign them as incorrect. We use a linear fit to find out outliers
    Parameters
    ----------
    bottom_shadow_pixels

    Returns
    -------

    """
    number_bar_positions, number_tubes = bottom_shadow_pixels.shape
    for tube_index in range(number_tubes):
        # Correct corner case 1
        jump_index = None
        y = bottom_shadow_pixels[:, tube_index]
        for i in range(number_bar_positions - 2, 0, -1):
            if abs(y[i] - y[i-1]) > 10 * max(1, abs(y[i+1] - y[i])):  # value/factor of 10 selected as threshold
                jump_index = i
                break
        if jump_index is not None:
            y[jump_index:] = INCORRECT_PIXEL_ASSIGNMENT
        # Correct corner case 2
        y = bottom_shadow_pixels[:, tube_index]
        x = bar_positions[y != INCORRECT_PIXEL_ASSIGNMENT]
        coefficients = np.polynomial.polynomial.polyfit(x, y[y != INCORRECT_PIXEL_ASSIGNMENT], 1)
        y_fitted = np.polynomial.polynomial.polyval(bar_positions, coefficients)
        residuals = np.abs(y - y_fitted)
        y[residuals > np.average(residuals) + 2 * np.std(residuals)] = INCORRECT_PIXEL_ASSIGNMENT


def calculate_apparent_tube_width(flood_input, component='detector1', load_barscan_calibration=True):
    r"""
    Determine the tube width most efficient for detecting neutrons. An effective tube (or pixel) diameter is
    determined for tubes in the front panel, and likewise for the tubes in the back panel.

    devs - Jose Borreguero <borreguerojm@ornl.gov>

    Parameters
    ----------
    flood_input: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Path to flood run, flood workspace name, or flood workspace object.
    component: str
        Name of the instrument component containing the detector array consisting of two parallel panels of tubes.
    load_barscan_calibration: bool
        Load pixel positions and heights from the pixel-calibrations database appropriate to ``input_workspace``. If
        py:obj:`False`, then the pixel positions and heigths will be those of ``input_workspace``.

    **Mantid algorithms used:**
        :ref:`DeleteWorkspaces <algm-DeleteWorkspaces-v1>`,
        :ref:`Integration <algm-Integration-v1>`,
        :ref:`MaskDetectors <algm-MaskDetectors-v1>`,
        :ref:`MaskDetectorsIf <algm-MaskDetectorsIf-v1>`,
        :ref:`ReplaceSpecialValues <algm-ReplaceSpecialValues-v1>`,

    Returns
    -------
    Returns
    -------
    dict
        Dictionary containing the following keys:
        - instrument, str, Standard name of the instrument.
        - component, str, name of the double detector array, usually "detector1".
        - run, int, run number of ``input_workspace``.
        - unit: str, the units for the tube widths. Set to 'mm' for mili-meters.
        - widths, list, A two-item list containing the apparent widths for the front and back tubes.
    """
    # Determine the type of input for the flood data
    if file_exists(flood_input):
        input_workspace = unique_workspace_dundername()
        Load(Filename=flood_input, OutputWorkspace=input_workspace)
    else:
        input_workspace = flood_input

    integrated_intensities = unique_workspace_dundername()
    Integration(InputWorkspace=input_workspace, OutputWorkspace=integrated_intensities)

    # Mask non-finite pixel_intensities (nan, inf). They can't be used in the calculation.
    #
    # Replace non-finite pixel_intensities with a value of -1
    ReplaceSpecialValues(InputWorkspace=integrated_intensities, OutputWorkspace=integrated_intensities,
                         NanValue=-1, NanError=-1, InfinityValue=-1, InfinityError=-1)
    # Mask detectors with negative pixel_intensities
    mask_workspace = unique_workspace_dundername()
    MaskDetectorsIf(InputWorkspace=integrated_intensities, Operator='Less', Value=0., OutputWorkspace=mask_workspace)
    MaskDetectors(Workspace=integrated_intensities, MaskedWorkspace=mask_workspace)

    # Update pixel positions and heights with the appropriate calibration, if so requested.
    if load_barscan_calibration is True:
        calibration = load_calibration(input_workspace, 'BARSCAN', component=component)
        calibration.apply(input_workspace)

    # Calculate the count density for each tube. Notice that if the whole tube is masked, then the associated
    # intensity is stored as nan.
    #
    # Sort the tubes according to the X-coordinate in decreasing value. This is the order when sitting on the
    # sample and iterating over the tubes "from left to right"
    collection = TubeCollection(integrated_intensities, 'detector1').sorted(view='decreasing X')
    detector_ids = list(itertools.chain.from_iterable(tube.detector_ids for tube in collection))
    count_densities = list()
    for tube in collection:
        weighted_intensities = tube.readY.ravel() / tube.pixel_heights
        d = np.mean(weighted_intensities[~tube.isMasked])
        count_densities.append(d)
    count_densities = np.array(count_densities)  # is convenient to cast densities into a numpy array data structure.

    # Determine the count densities per panel and for the whole detector array.
    # We must be careful to pick only tubes with finite densities (avoid 'nan')
    average_count_density = np.mean(count_densities[np.isfinite(count_densities)])
    front_count_density = np.mean(count_densities[::2][np.isfinite(count_densities[::2])])  # front tubes, even indexes
    back_count_density = np.mean(count_densities[1::2][np.isfinite(count_densities[1::2])])  # back tubes, odd indexes

    # Determine the front and back pixel widths
    nominal_width = collection[0][0].width  # width of the first pixel in the first tube
    front_width = (front_count_density / average_count_density) * nominal_width
    back_width = (back_count_density / average_count_density) * nominal_width

    # Generate a list of pixel widths. It is assumed that front tubes have an even tube index
    widths = list()
    for tube_index, tube in enumerate(collection):
        pixel_width = front_width if tube_index % 2 == 0 else back_width
        widths.extend([pixel_width] * len(tube))

    DeleteWorkspaces(integrated_intensities, mask_workspace)

    metadata = dict(caltype=CalType.TUBEWIDTH.name,
                    instrument=instrument_standard_name(input_workspace),
                    component=component,
                    daystamp=day_stamp(input_workspace),
                    runnumbers=[SampleLogs(input_workspace).single_value('run_number'), ])
    return Table(metadata, detector_ids, widths=widths)
