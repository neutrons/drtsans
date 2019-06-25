"""
    Standalone mock api to be used as a prototype example for
    instrument scientists. Those function should be put in the real code
    and replaced as we implement newer versions.
"""
from __future__ import (absolute_import, division, print_function)
import logging
import mantid.simpleapi as api
import os


def find_beam_center(filename):
    """
        Find beam center

        :param ws: scattering workspace
        :param filename: file path for direct beam data
    """
    x, y, _ = api.SANSBeamFinder(Filename=filename)

    return x, y


def load_events(filename, beam_center_x=None, beam_center_y=None,
                tof_min_cut=0, tof_max_cut=0,
                sample_offset=0):
    """
        Load data should load, move the detector in the right place,
        and correct the TOF.
        For EQSANS, this is currently EQSANSLoad

        @param tof_min_cut: amount of time to cut at the beginning
                            of distribution [musec]
        @param tof_max_cut: amount of time to cut at the end
                            of distribution [musec]
        @param sample_offset: Offset to be applied to the sample position [mm]
    """
    _, output_ws = os.path.split(filename)
    ws, _, _ = api.EQSANSLoad(Filename=filename,
                              BeamCenterX=beam_center_x,
                              BeamCenterY=beam_center_y,
                              UseConfigTOFCuts=False,
                              SampleOffset=sample_offset,
                              DetectorOffset=0,
                              LowTOFCut=tof_min_cut, HighTOFCut=tof_max_cut,
                              SkipTOFCorrection=False, WavelengthStep=0.1,
                              UseConfigMask=False, OutputWorkspace=output_ws)
    return ws


def apply_transmission(workspace, beam_center_x=None, beam_center_y=None):
    """
        Computes transmission and error and applies it.
        :param beam_center_x: beam center in x for transmission run (in pixel)
        :param beam_center_y: beam center in y for transmission run (in pixel)

        Note: comes from this algorithm
            # p=property_manager.getProperty("TransmissionAlgorithm")
    """
    logging.warning("Not yet implemented")
    return workspace


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


def absolute_scale(ws):
    """
        Normalize to absolute scale

        :param ws: scattering workspace

        Note: comes from this algorithm
           self._simple_execution("AbsoluteScaleAlgorithm", output_ws)
    """
    logging.warning("Not yet implemented")
    return ws


def geometry_correction(ws):
    """
        Geometry correction

        Note: comes from the following algorithm
            self._simple_execution("GeometryAlgorithm", output_ws)
    """
    logging.warning("Not yet implemented")
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
        return outputs[0], outputs[1]
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
                 beam_profile=''):
    """
        Prepare data for reduction
        :param ws: scattering workspace
        :param normalize_to_monitor: if False, normalize to proton charge
        :param beam_profile: if supplied, divide by the beam profile
    """
    # Dark current subtraction
    # self._simple_execution("DarkCurrentAlgorithm", workspace)

    # Normalize
    normalize_to_beam = os.path.isfile(beam_profile)
    ws, _ = api.EQSANSNormalise(InputWorkspace=workspace,
                                NormaliseToBeam=normalize_to_beam,
                                BeamSpectrumFile=beam_profile,
                                NormaliseToMonitor=normalize_to_monitor)

    # Mask
    api.SANSMask(Workspace=ws, MaskedDetectorList=None,
                 MaskedEdges=[0, 0, 0, 0], MaskedSide=None)

    # Solid angle correction
    ws, _ = api.SANSSolidAngleCorrection(InputWorkspace=ws,
                                         DetectorTubes=True,
                                         DetectorWing=False,
                                         OutputWorkspace=ws)

    # Sensitivity correction
    # p=property_manager.getProperty("SensitivityBeamCenterAlgorithm")
    # p=property_manager.getProperty("SensitivityAlgorithm")
    return ws
