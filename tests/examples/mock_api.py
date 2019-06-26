# flake8: noqa
"""
    Standalone mock api to be used as a prototype example for
    instrument scientists. Those function should be put in the real code
    and replaced as we implement newer versions.
"""
from __future__ import (absolute_import, division, print_function)
import mantid.simpleapi as api
import os


def find_beam_center(filename):
    """
        Find beam center

        :param filename: file path for direct beam data
    """
    x, y, _ = api.SANSBeamFinder(Filename=filename)
    return x, y


def load_events(filename, beam_center_x=None, beam_center_y=None,
                tof_min_cut=None, tof_max_cut=None, use_config_tof_cuts=False,
                sample_offset=0, use_config=False, use_config_mask=False,
                correct_for_flight_path=False):
    """
        Load data should load, move the detector in the right place,
        and correct the TOF.
        For EQSANS, this is currently EQSANSLoad

        :param filename: data file to load
        :param beam_center_x: beam center in x for transmission run (in pixel)
        :param beam_center_y: beam center in y for transmission run (in pixel)
        :param tof_min_cut: amount of time to cut at the beginning
                            of distribution [musec]
        :param tof_max_cut: amount of time to cut at the end
                            of distribution [musec]
        :param use_config_tof_cuts: if true, apply TOF cuts according to config
        :param use_config: if true, use EQSANS config file
        :param use_config_mask: if true, apply mask according to config
        :param correct_for_flight_path: if true, correct for flight to detector
        :param sample_offset: Offset to be applied to the sample position [mm]
    """
    _, output_ws = os.path.split(filename)
    ws, _ = api.EQSANSLoad(Filename=filename,
                           BeamCenterX=beam_center_x,
                           BeamCenterY=beam_center_y,
                           UseConfigTOFCuts=use_config_tof_cuts,
                           SampleOffset=sample_offset,
                           DetectorOffset=0,
                           LoadMonitors=False,
                           PreserveEvents=False,
                           UseConfig=use_config,
                           CorrectForFlightPath=correct_for_flight_path,
                           LowTOFCut=tof_min_cut, HighTOFCut=tof_max_cut,
                           SkipTOFCorrection=False, WavelengthStep=0.1,
                           UseConfigMask=use_config_mask,
                           OutputWorkspace=output_ws)
    return ws

def apply_transmission(workspace, transmission_value=1.0,
                       transmission_error=0.0,
                       theta_dependent=False):
    """
        Apply a simple transmission
        :param workspace: workspace to apply the correction to
        :param transmission_value: transmission value
        :param transmission_error: transmission error
    """
    ws = api.ApplyTransmissionCorrection(InputWorkspace=workspace,
                                         TransmissionValue=transmission_value,
                                         TransmissionError=transmission_error,
                                         ThetaDependent=theta_dependent,
                                         OutputWorkspace=str(workspace))
    return ws

def compute_and_apply_transmission(workspace,
                                   sample_data_file=None, empty_data_file=None,
                                   beam_center_x=None, beam_center_y=None,
                                   dark_data_file=None):
    """
        Computes transmission and error and applies it.
        :param workspace: workspace to apply the correction to
        :param sample_data_file: sample direct beam file
        :param empty_data_file: empty direct beam file
        :param beam_center_x: beam center in x for transmission run (in pixel)
        :param beam_center_y: beam center in y for transmission run (in pixel)
        :param dark_data_file: dark current data file

        Note: This should be broken down in two algorithms
    """
    ws, _, _, _ = api.EQSANSDirectBeamTransmission(InputWorkspace=workspace,
                                                   SampleDataFilename=sample_data_file,
                                                   EmptyDataFilename=empty_data_file,
                                                   DarkCurrentFilename=dark_data_file,
                                                   UseSampleDarkCurrent=False,
                                                   BeamRadius=5,
                                                   ThetaDependent=True,
                                                   FitFramesTogether=False,
                                                   BeamCenterX=beam_center_x,
                                                   BeamCenterY=beam_center_y,
                                                   TransmissionWorkspace='tr_ws',
                                                   RawTransmissionWorkspace='raw_tr',
                                                   OutputWorkspace=str(workspace))
    return ws


def subtract_dark_current(workspace, dark_data_file):
    """
        Subtract dark current
        Note: in the future we may want to do something simpler and
        pass a workspace

        :param workspace: workspace to apply the correction to
        :param dark_data_file: dark current data file
    """
    outputs = api.EQSANSDarkCurrentSubtraction(InputWorkspace=workspace,
                                               Filename=dark_data_file,
                                               OutputDarkCurrentWorkspace='dark_current',
                                               OutputWorkspace=str(workspace))
    return outputs[0]


def subtract_background(ws, ws_bck):
    """
        Subtract background.
        The reason this is a function instead of just the minus
        operator is that we need a rebin in between.

        :param ws: scattering workspace
        :param ws: background workspace
    """
    rebinned_ws = '%s_rebin' % ws_bck
    api.RebinToWorkspace(WorkspaceToRebin=ws_bck,
                         WorkspaceToMatch=ws,
                         OutputWorkspace=rebinned_ws,
                         PreserveEvents=True)
    ws = api.Minus(LHSWorkspace=ws,
                   RHSWorkspace=rebinned_ws,
                   OutputWorkspace=ws)
    return ws


def iq(ws, number_of_bins=100, log_binning=False, sample_aperture=10.0):
    """
        Compute I(q)
        :param ws: scattering workspace
        :param number_of_bins int: number of bins in Q
        :param log_binning bool: if True, use log binning in Q
        :param sample_aperture float: sample aperture in mm for resolution
    """
    outputs = api.EQSANSAzimuthalAverage1D(InputWorkspace=ws,
                                           NumberOfBins=number_of_bins,
                                           LogBinning=log_binning,
                                           IndependentBinning=True,
                                           ScaleResults=True,
                                           ComputeResolution=True,
                                           SampleApertureDiameter=sample_aperture,  # noqa: E501
                                           OutputWorkspace="%s_iq" % ws)

    # For frame-skipping mode, we have an additional workspace do deal with
    if len(outputs) == 3:
        return outputs[0], outputs[2]
    else:
        return outputs[0]


def iqxqy(ws, number_of_bins=100, log_binning=False):
    """
        Compute I(qx,qy)
        :param ws: scattering workspace
        :param number_of_bins int: number of bins in Q
        :param log_binning bool: if True, use log binning in Q
    """
    iq_ws, _, _ = api.EQSANSQ2D(InputWorkspace=ws,
                                NumberOfBins=number_of_bins,
                                IQxQyLogBinning=log_binning,
                                OutputWorkspace="%s_iqxy" % ws)
    return iq_ws


def prepare_data(workspace, normalize_to_monitor=False,
                 beam_center_x=None, beam_center_y=None,
                 beam_profile='', tubes=True,
                 dark_data_file=None, masked_detector_list=None,
                 sensitivity_file=None):
    """
        Prepare data for reduction
        :param ws: scattering workspace
        :param normalize_to_monitor: if False, normalize to proton charge
        :param beam_center_x: beam center in x for transmission run (in pixel)
        :param beam_center_y: beam center in y for transmission run (in pixel)
        :param beam_profile: if supplied, divide by the beam profile
        :param dark_data_file: dark current data file
        :param masked_detector_list: additional pixels to mask
        :param sensitivity_file: processed sensitivity file
        :param tubes: whether of not we have detector tubes
    """
    # Dark current subtraction
    if dark_data_file is not None:
        workspace = subtract_dark_current(workspace, dark_data_file)

    # Normalize
    normalize_to_beam = os.path.isfile(beam_profile)
    ws, _ = api.EQSANSNormalise(InputWorkspace=workspace,
                                NormaliseToBeam=normalize_to_beam,
                                BeamSpectrumFile=beam_profile,
                                NormaliseToMonitor=normalize_to_monitor,
                                OutputWorkspace=str(workspace))

    # Mask
    api.SANSMask(Workspace=ws, MaskedDetectorList=masked_detector_list,
                 MaskedEdges=[0, 0, 0, 0], MaskedSide=None)

    # Solid angle correction
    ws, _ = api.SANSSolidAngleCorrection(InputWorkspace=ws,
                                         DetectorTubes=tubes,
                                         DetectorWing=False,
                                         OutputWorkspace=str(ws))

    # Sensitivity correction
    if sensitivity_file is not None:
        ws, _, _ = api.SANSSensitivityCorrection(InputWorkspace=ws,
                                                 Filename=sensitivity_file,
                                                 DarkCurrentFile=dark_data_file,
                                                 MinEfficiency=0.5, MaxEfficiency=1.5,
                                                 BeamCenterX=beam_center_x,
                                                 BeamCenterY=beam_center_y,
                                                 OutputSensitivityWorkspace='sensitivity',
                                                 OutputWorkspace=str(ws))

    return ws
