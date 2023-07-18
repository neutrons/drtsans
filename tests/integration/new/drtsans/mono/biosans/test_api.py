# local imports
from drtsans.mono.biosans.beam_finder import find_beam_center
from drtsans.mono.biosans.api import load_all_files, process_single_configuration, reduce_single_configuration
from drtsans.mono.load import transform_to_wavelength
from drtsans.mono.transmission import calculate_transmission
from drtsans.redparms import reduction_parameters
from drtsans.settings import unique_workspace_dundername

# third party imports
from mantid.api import AnalysisDataService
from mantid.dataobjects import Workspace2D
from mantid.simpleapi import (
    mtd,
    CloneWorkspace,
    DeleteWorkspace,
    DeleteWorkspaces,
    ExtractMask,
    InvertMask,
    LoadNexusProcessed,
    MaskAngle,
    MaskBTP,
)
import numpy as np
from numpy.testing import assert_allclose

import pytest

# standard imports
from math import isclose
import os
from os.path import join as pjoin
import threading
from typing import Union
from unittest.mock import patch as mock_patch


def _mock_LoadEventNexus(*args, **kwargs):
    # Substitute LoadEventNexus with LoadNexusProcessed because our synthetic files were created with SaveNexus
    return LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])


# need to hardcode these worspace removals, is there a way to derive them?
# is there a safe way to just ping for a list of workspaces generated by this user?
workspaces = [
    "_bkgd_trans",
    "_empty",
    "_filter",
    "_info",
    "_load_tmp",
    "_load_tmp_monitors",
    "_main_sensitivity",
    "_sample_trans",
    "_wing_sensitivity",
    "barscan_BIOSANS_detector1_20210302",
    "barscan_BIOSANS_wing_detector_20210302",
    "processed_data_main_0",
    "processed_data_main_1",
    "processed_data_main_10",
    "processed_data_main_2",
    "processed_data_main_3",
    "processed_data_main_4",
    "processed_data_main_5",
    "processed_data_main_6",
    "processed_data_main_7",
    "processed_data_main_8",
    "processed_data_main_9",
    "processed_data_wing_0",
    "processed_data_wing_1",
    "processed_data_wing_10",
    "processed_data_wing_2",
    "processed_data_wing_3",
    "processed_data_wing_4",
    "processed_data_wing_5",
    "processed_data_wing_6",
    "processed_data_wing_7",
    "processed_data_wing_8",
    "processed_data_wing_9",
    "TOFCorrectWS",
    "tubewidth_BIOSANS_detector1_20200613",
    "tubewidth_BIOSANS_detector1_20210304",
    "tubewidth_BIOSANS_wing_detector_20200613",
    "tubewidth_BIOSANS_wing_detector_20210304",
]


def test_process_single_configuration(biosans_synthetic_dataset, clean_workspace):
    r"""
    Apply process_single_configuration() to the synthetic sample data. This will correct the sample data, among other
    things, by subtracting the background (which is also a synthetic data).
    We check that the background data is correctly subtracted by looking in the sample data for the average intensity
    of the pixels in two regions of the midrange detector. Were not for the background, the average intensity in these
    two regions would be the same.
    The first region is a ring sector within the Midrange detector defined by lower and upper values in scattering
    angle twotheta (5.70 and 6.20 degrees, respectively).
    The second region is also a ring sector within the Midrange detector defined by lower and upper values in
    scattering angle twotheta (8.64 and 8.99 degrees, respectively).
    The background data affects only the first region. As a consequence, the average intensity in the first region
    is about 2.5 times that of the second region.
    After applying process_single_configuration() to the sample data (which includes the background subtraction),
    the average intensity in the first region is about 1.0 times that of the second region, as expected.
    """

    def _ring_maskworkspace(
        input_workspace: Union[str, Workspace2D], twotheta_begin: float, twotheta_end: float
    ) -> str:
        r"""
        Create a MaskWorkspace where everything is masked except for a ring sector within the Midrange detector. A
        ring is defined by a lower and upper boundary in scattering angle twotheta.

        Parameters
        ----------
        input_workspace
             Workspace to be used as template for the output MaskWorkspace
        twotheta_begin
            lower boundary of the ring sector in degrees
        twotheta_end
            upper boundary of the ring sector in degrees
        Returns
        -------
        Name of the MaskWorkspace
        """
        ring_mask = unique_workspace_dundername()
        CloneWorkspace(InputWorkspace=input_workspace, OutputWorkspace=ring_mask)
        MaskAngle(Workspace=ring_mask, MinAngle=twotheta_begin, MaxAngle=twotheta_end, Angle="TwoTheta")
        ExtractMask(InputWorkspace=ring_mask, OutputWorkspace=ring_mask)
        InvertMask(InputWorkspace=ring_mask, OutputWorkspace=ring_mask)
        MaskBTP(Workspace=ring_mask, Components="detector1,wing_detector")
        clean_workspace(ring_mask)  # mark for deletion
        return ring_mask

    def _intensity_ratio(sample_workspace: Union[str, Workspace2D], first_mask: str, second_mask: str) -> float:
        r"""
        Given two mask workspaces, calculate the ratio of the average intensity of the unmasked pixels in the first
        mask to the average intensity of the unmasked pixels in the second mask.

        Parameters
        ----------
        sample_workspace
        first_mask
        second_mask
        """
        MASK_VALUE = 1
        intensities_all = mtd[str(sample_workspace)].extractY()
        average_intensities = list()  # average intensity of a pixel in the first and second ring sector
        for mask_workspace in (first_mask, second_mask):
            pixels_unmasked = np.where(mtd[mask_workspace].extractY() < MASK_VALUE)[0]
            average_intensities.append(np.average(intensities_all[pixels_unmasked]))
        return average_intensities[0] / average_intensities[1]

    def _load_synthetic_run(keyword: str) -> Workspace2D:
        r"""Load a synthetic run from the biosans_synthetic_dataset fixture and histogram to wavelength"""
        workspace = LoadNexusProcessed(
            Filename=pjoin(
                biosans_synthetic_dataset["data_dir"], f"CG3_{biosans_synthetic_dataset['runs'][keyword]}.nxs.h5"
            ),
            OutputWorkspace=unique_workspace_dundername(),
        )
        workspace = transform_to_wavelength(workspace)
        clean_workspace(workspace)  # mark it for deletion
        return workspace

    # Find the beam center
    ws_beam_center = _load_synthetic_run("beam_center")
    center_x, center_y, center_y_wing, center_y_midrange, _ = find_beam_center(
        ws_beam_center, centering_options={"IntegrationRadius": 0.03}
    )

    workspace_sample = _load_synthetic_run("sample")
    first_ring_mask = _ring_maskworkspace(workspace_sample, 5.7, 6.2)
    second_ring_mask = _ring_maskworkspace(workspace_sample, 8.64, 8.99)
    # ratio of average pixel intensity in the first ring sector to that average intensity in the second ring sector
    ratio = _intensity_ratio(workspace_sample, first_ring_mask, second_ring_mask)
    assert_allclose(ratio, 2.56, 0.1)

    prefix = "ewjsfvshisn"  # random string prepended to all workspaces produced by the function
    process_single_configuration(
        workspace_sample,
        sample_trans_ws=_load_synthetic_run("sample_transmission"),
        bkg_ws_raw=_load_synthetic_run("background"),
        bkg_trans_ws=_load_synthetic_run("background_transmission"),
        center_x=center_x,
        center_y=center_y,
        center_y_wing=center_y_wing,
        center_y_midrange=center_y_midrange,
        dark_current=_load_synthetic_run("dark_current"),
        flux_method="time",
        output_prefix=prefix,
    )

    # verify both ring sectors have the same average intensity per pixel
    ratio = _intensity_ratio(prefix + "_sample", first_ring_mask, second_ring_mask)
    assert_allclose(ratio, 1.0, 0.1)

    DeleteWorkspaces([prefix + "_" + suffix for suffix in ("sample", "background")])


@pytest.mark.skipif(
    not os.path.exists("/HFIR/HB2B/shared/autoreduce/"),
    reason="Skip test on build server",
)
def test_reduce_single_configuration_slice_transmission_false(temp_directory):
    reduction_input = {
        "schemaStamp": "2020-04-15T21:09:52.745905",
        "instrumentName": "BIOSANS",
        "iptsNumber": "24666",
        "dataDirectories": None,
        "sample": {
            "runNumber": "8375",
            "thickness": 0.2,
            "transmission": {"runNumber": "8379", "value": None},
        },
        "background": {
            "runNumber": "8374",
            "transmission": {"runNumber": "8378", "value": None},
        },
        "emptyTransmission": {"runNumber": "8381", "value": None},
        "beamCenter": {"runNumber": "8381"},
        "outputFileName": "r8375_AgBeh_15m18Aqa",
        "configuration": {
            "wavelength": None,
            "wavelengthSpread": None,
            "useTimeSlice": False,
            "useTimeSliceTransmission": False,
            "timeSliceInterval": 60.0,
            "useLogSlice": False,
            "logSliceName": None,
            "logSliceInterval": None,
            "sampleOffset": None,
            "sampleApertureSize": 12.0,
            "sampleDetectorDistance": None,
            "sampleToSi": None,
            "sourceApertureDiameter": None,
            "usePixelCalibration": True,
            "maskFileName": None,
            "useDefaultMask": True,
            "defaultMask": [
                {"Pixel": "1-18,239-256"},
                {"Bank": "18-24,42-48"},
                {"Bank": "49", "Tube": "1"},
                {"Bank": "88", "Tube": "4"},
            ],
            "useMaskBackTubes": False,
            "darkMainFileName": "CG3_8331.nxs.h5",
            "darkWingFileName": "CG3_8331.nxs.h5",
            "darkMidrangeFileName": "CG3_8331.nxs.h5",
            "normalization": "Monitor",
            "normalizationResortToTime": False,
            "sensitivityMainFileName": "/HFIR/CG3/shared/Cycle490/Sens_f8367m7p0_bsSVP.nxs",
            "sensitivityWingFileName": "/HFIR/CG3/shared/Cycle490/Sens_f8369w1p4_bsSVP.nxs",
            # sensitivity file for mid-range detector should be added after
            # the installation of the new mid-range detector.
            "useSolidAngleCorrection": True,
            "blockedBeamRunNumber": None,
            "useThetaDepTransCorrection": True,
            "DBScalingBeamRadius": 40.0,
            "mmRadiusForTransmission": None,
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": 2.286e-09,
            "numMainQxQyBins": 100,
            "numWingQxQyBins": 100,
            "1DQbinType": "scalar",
            "QbinType": "log",
            "numMainQBins": None,
            "numWingQBins": None,
            "LogQBinsPerDecadeMain": 25,
            "LogQBinsPerDecadeWing": 25,
            "useLogQBinsDecadeCenter": False,
            "useLogQBinsEvenDecade": False,
            "WedgeMinAngles": None,
            "WedgeMaxAngles": None,
            "autoWedgeQmin": 0.003,
            "autoWedgeQmax": 0.04,
            "autoWedgeQdelta": 0.01,
            "autoWedgeAzimuthalDelta": 1.0,
            "autoWedgePeakWidth": 0.5,
            "autoWedgeBackgroundWidth": 1.0,
            "autoWedgeSignalToNoiseMin": 2.0,
            "AnnularAngleBin": 1.0,
            "useErrorWeighting": False,
            "smearingPixelSizeX": None,
            "smearingPixelSizeY": None,
            "useSubpixels": False,
            "subpixelsX": None,
            "subpixelsY": None,
            "QminMain": 0.0009,
            "QmaxMain": 0.016,
            "QminWing": 0.009,
            "QmaxWing": 0.3,
            "overlapStitchQmin": [0.0105],
            "overlapStitchQmax": [0.0145],
            "wedge1QminMain": 0.003,
            "wedge1QmaxMain": 0.0425,
            "wedge1QminWing": 0.02,
            "wedge1QmaxWing": 0.45,
            "wedge1overlapStitchQmin": 0.025,
            "wedge1overlapStitchQmax": 0.04,
            "wedge2QminMain": 0.003,
            "wedge2QmaxMain": 0.0425,
            "wedge2QminWing": 0.03,
            "wedge2QmaxWing": 0.45,
            "wedge2overlapStitchQmin": 0.03,
            "wedge2overlapStitchQmax": 0.04,
            "wedges": None,
            "symmetric_wedges": True,
        },
        "logslice_data": {},
    }
    reduction_input["configuration"]["outputDir"] = temp_directory(prefix="trans_slice_false")

    prefix = "sans-backend-test" + str(threading.get_ident()) + "_"
    loaded = load_all_files(reduction_input, prefix)
    _ = reduce_single_configuration(loaded, reduction_input)
    # just need a couple components from reduce
    # but the whole thing needs to be run then a few components pulled
    transmission = calculate_transmission(
        mtd["_sample_trans"],  # pull relevant transmission
        mtd["_empty"],  # pull relevant
    )
    transmission_val = transmission.extractY()[0][0]
    clean_all_ws(prefix)
    assert isclose(transmission_val, 0.5734218305525239)  # provided by s6v
    del _


def clean_all_ws(prefix):
    for workspace in workspaces:
        remove_ws(workspace)
    for object_name in mtd.getObjectNames():
        if object_name.startswith(prefix):
            remove_ws(object_name)


def remove_ws(workspace):
    if AnalysisDataService.doesExist(workspace):
        DeleteWorkspace(workspace)


def test_reduce_single_configuration_slice_transmission_true(temp_directory):
    reduction_input = {
        "schemaStamp": "2020-04-15T21:09:52.745905",
        "instrumentName": "BIOSANS",
        "iptsNumber": "24666",
        "dataDirectories": None,
        "sample": {
            "runNumber": "8361",
            "thickness": 0.1,
            "transmission": {"runNumber": "8361", "value": None},
        },
        "background": {
            "runNumber": "8359",
            "transmission": {"runNumber": "8359", "value": None},
        },
        "emptyTransmission": {"runNumber": "8364", "value": None},
        "beamCenter": {"runNumber": "8373"},
        "IntegrationRadius": None,
        "outputFileName": "r8361_PorB3_15m",
        "configuration": {
            "wavelength": None,
            "wavelengthSpread": None,
            "useTimeSlice": True,
            "useTimeSliceTransmission": True,
            "timeSliceInterval": 60.0,
            "useLogSlice": False,
            "logSliceName": None,
            "logSliceInterval": None,
            "sampleOffset": None,
            "sampleApertureSize": 14.0,
            "sampleDetectorDistance": None,
            "sampleToSi": None,
            "sourceApertureDiameter": None,
            "usePixelCalibration": True,
            "maskFileName": None,
            "useDefaultMask": True,
            "defaultMask": [
                {"Pixel": "1-18,239-256"},
                {"Bank": "18-24,42-48"},
                {"Bank": "49", "Tube": "1"},
                {"Bank": "88", "Tube": "4"},
            ],
            "useMaskBackTubes": False,
            "darkMainFileName": "CG3_8331.nxs.h5",
            "darkWingFileName": "CG3_8331.nxs.h5",
            "darkMidrangeFileName": "CG3_8331.nxs.h5",
            "normalization": "Monitor",
            "normalizationResortToTime": False,
            "sensitivityMainFileName": "/HFIR/CG3/shared/Cycle490/Sens_f8367m7p0_bsSVP.nxs",
            "sensitivityWingFileName": "/HFIR/CG3/shared/Cycle490/Sens_f8369w1p4_bsSVP.nxs",
            # sensitivity file for mid-range detector should be added after
            # the installation of the new mid-range detector.
            "useSolidAngleCorrection": True,
            "blockedBeamRunNumber": None,
            "useThetaDepTransCorrection": True,
            "DBScalingBeamRadius": 40.0,
            "mmRadiusForTransmission": None,
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": 2.094e-10,
            "numMainQxQyBins": 100,
            "numWingQxQyBins": 100,
            "1DQbinType": "scalar",
            "QbinType": "log",
            "numMainQBins": None,
            "numWingQBins": None,
            "LogQBinsPerDecadeMain": 25,
            "LogQBinsPerDecadeWing": 25,
            "useLogQBinsDecadeCenter": False,
            "useLogQBinsEvenDecade": False,
            "WedgeMinAngles": None,
            "WedgeMaxAngles": None,
            "autoWedgeQmin": 0.003,
            "autoWedgeQmax": 0.04,
            "autoWedgeQdelta": 0.01,
            "autoWedgeAzimuthalDelta": 1.0,
            "autoWedgePeakWidth": 0.5,
            "autoWedgeBackgroundWidth": 1.0,
            "autoWedgeSignalToNoiseMin": 2.0,
            "AnnularAngleBin": 1.0,
            "useErrorWeighting": False,
            "smearingPixelSizeX": None,
            "smearingPixelSizeY": None,
            "useSubpixels": False,
            "subpixelsX": None,
            "subpixelsY": None,
            "QminMain": 0.003,
            "QmaxMain": 0.045,
            "QminWing": 0.03,
            "QmaxWing": 0.9,
            "overlapStitchQmin": [0.0325],
            "overlapStitchQmax": [0.0425],
            "wedge1QminMain": 0.003,
            "wedge1QmaxMain": 0.0425,
            "wedge1QminWing": 0.02,
            "wedge1QmaxWing": 0.45,
            "wedge1overlapStitchQmin": 0.025,
            "wedge1overlapStitchQmax": 0.04,
            "wedge2QminMain": 0.003,
            "wedge2QmaxMain": 0.0425,
            "wedge2QminWing": 0.03,
            "wedge2QmaxWing": 0.45,
            "wedge2overlapStitchQmin": 0.03,
            "wedge2overlapStitchQmax": 0.04,
            "wedges": None,
            "symmetric_wedges": True,
        },
        "logslice_data": {},
    }
    reduction_input["configuration"]["outputDir"] = temp_directory(prefix="trans_slice_true")
    prefix = "sans-backend-test" + str(threading.get_ident()) + "_"
    loaded = load_all_files(reduction_input, prefix)
    _ = reduce_single_configuration(loaded, reduction_input)
    # just need a couple components from reduce
    # but the whole thing needs to be run then a few components pulled
    transmission = calculate_transmission(
        mtd["_sample_trans"],  # pull relevant transmission
        mtd["_empty"],  # pull relevant
    )

    transmission_val = transmission.extractY()[0][0]
    clean_all_ws(prefix)
    assert isclose(transmission_val, 0.7526460467895154)  # from above config using older workflow
    del _


@pytest.mark.skip(reason="Skip until the mantid nightly conda package is updated, with timestamp > 2030-07-12")
@mock_patch("drtsans.load.LoadEventNexus", new=_mock_LoadEventNexus)
@mock_patch("drtsans.load.__monitor_counts")
def test_reduce_single_configuration_synthetic_dataset(mock_monitor_counts, biosans_synthetic_dataset, temp_directory):
    data = biosans_synthetic_dataset
    mock_monitor_counts.return_value = biosans_synthetic_dataset["monitor_counts"]
    reduction_input = {
        "schemaStamp": "2020-04-15T21:09:52.745905",
        "instrumentName": "CG3",
        "iptsNumber": "00000",
        "dataDirectories": f"{data['data_dir']}",
        "sample": {
            "runNumber": "92310",
            "thickness": 0.2,
            "transmission": {"runNumber": "92330", "value": None},
        },
        "background": {
            "runNumber": "92320",
            "transmission": {"runNumber": "92330", "value": None},
        },
        "emptyTransmission": {"runNumber": "92300", "value": None},
        "beamCenter": {
            "runNumber": "92300",
            "method": "center_of_mass",
            "com_centering_options": {"IntegrationRadius": 0.07},
        },
        "outputFileName": "synthethic",
        "configuration": {
            "outputDir": temp_directory(prefix="sythetic_experiment"),
            "timeSliceInterval": 60.0,
            "sampleApertureSize": 12.0,
            "usePixelCalibration": True,
            "useDefaultMask": True,
            "defaultMask": [
                {"Pixel": "1-18,239-256"},
            ],
            "darkMainFileName": "CG3_92340.nxs.h5",
            "darkWingFileName": "CG3_92340.nxs.h5",
            "darkMidrangeFileName": "CG3_92340.nxs.h5",
            "sensitivityMainFileName": os.path.join(data["data_dir"], "sensitivity_detector1.nxs"),
            "sensitivityWingFileName": os.path.join(data["data_dir"], "sensitivity_wing_detector.nxs"),
            "sensitivityMidrangeFileName": os.path.join(data["data_dir"], "sensitivity_midrange_detector.nxs"),
            "useThetaDepTransCorrection": True,
            "DBScalingBeamRadius": 40.0,
            "StandardAbsoluteScale": 2.286e-09,
            "numMainQxQyBins": 100,
            "numWingQxQyBins": 100,
            "1DQbinType": "scalar",
            "QbinType": "log",
            "LogQBinsPerDecadeMain": 25,
            "LogQBinsPerDecadeWing": 25,
            "WedgeMinAngles": None,
            "WedgeMaxAngles": None,
            "autoWedgeQmin": 0.003,
            "autoWedgeQmax": 0.04,
            "autoWedgeQdelta": 0.01,
            "autoWedgeAzimuthalDelta": 1.0,
            "autoWedgePeakWidth": 0.5,
            "autoWedgeBackgroundWidth": 1.0,
            "autoWedgeSignalToNoiseMin": 2.0,
            "AnnularAngleBin": 1.0,
            "QminMain": 0.0009,
            "QmaxMain": 0.016,
            "QminWing": 0.009,
            "QmaxWing": 0.3,
            "overlapStitchQmin": [0.0105],
            "overlapStitchQmax": [0.0145],
            "wedge1QminMain": 0.003,
            "wedge1QmaxMain": 0.0425,
            "wedge1QminWing": 0.02,
            "wedge1QmaxWing": 0.45,
            "wedge1overlapStitchQmin": 0.025,
            "wedge1overlapStitchQmax": 0.04,
            "wedge2QminMain": 0.003,
            "wedge2QmaxMain": 0.0425,
            "wedge2QminWing": 0.03,
            "wedge2QmaxWing": 0.45,
            "wedge2overlapStitchQmin": 0.03,
            "wedge2overlapStitchQmax": 0.04,
            "wedges": None,
            "symmetric_wedges": True,
        },
        "logslice_data": {},
    }
    reduction_input = reduction_parameters(parameters_particular=reduction_input, validate=False)
    prefix = "sans-backend-test" + str(threading.get_ident()) + "_"
    loaded = load_all_files(reduction_input, prefix, path=data["data_dir"])
    reduce_single_configuration(loaded, reduction_input)


if __name__ == "__main__":
    pytest.main(__file__)
