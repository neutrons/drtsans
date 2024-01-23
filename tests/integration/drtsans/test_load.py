import os
import json
import pytest

from drtsans.load import load_events
from mantid.simpleapi import DeleteWorkspace, mtd


class TestLoadEvents:
    def test_pixel_calibration(self, datarepo_dir, temp_directory):
        """Check the pixel calibration is applied to a workspace upon loading"""

        # create pixel_calibration.json with tablefile pointing to data repo
        tmp_dir = temp_directory()
        pixel_calibration = [
            {
                "caltype": "BARSCAN",
                "instrument": "GPSANS",
                "component": "detector1",
                "daystamp": 20200103,
                "runnumbers": [7465],
                "tablefile": os.path.join(datarepo_dir.gpsans, "barscan_GPSANS_detector1_20200103.nxs"),
            },
            {
                "caltype": "TUBEWIDTH",
                "instrument": "GPSANS",
                "component": "detector1",
                "daystamp": 20200130,
                "runnumbers": [8143],
                "tablefile": os.path.join(datarepo_dir.gpsans, "tubewidth_GPSANS_detector1_20200130.nxs"),
            },
        ]
        pixel_calibration_filename = os.path.join(tmp_dir, "pixel_calibration.json")
        with open(pixel_calibration_filename, "w") as f:
            json.dump(pixel_calibration, f)

        file_name = os.path.join(datarepo_dir.gpsans, "CG2_8148.nxs.h5")
        workspace = load_events(
            file_name,
            pixel_calibration=pixel_calibration_filename,
            output_workspace=mtd.unique_hidden_name(),
        )
        component_info = workspace.componentInfo()
        component_info_index = 42  # identifies some detector pixel

        # Assert pixel width is calibrated
        nominal_pixel_width = component_info.shape(component_info_index).getBoundingBox().width().X()
        assert nominal_pixel_width == pytest.approx(0.00804, abs=1.0e-5)  # uncalibrated width
        pixel_width = component_info.scaleFactor(42).X() * nominal_pixel_width
        assert pixel_width == pytest.approx(0.00968, abs=1.0e-5)  # calibrated width

        # Assert pixel height is calibrated
        nominal_pixel_height = component_info.shape(component_info_index).getBoundingBox().width().Y()
        assert nominal_pixel_height == pytest.approx(0.00409, abs=1.0e-5)  # uncalibrated width
        pixel_height = component_info.scaleFactor(42).X() * nominal_pixel_height
        assert pixel_height == pytest.approx(0.00492, abs=1.0e-5)  # calibrated width

        # NOTE:
        # These workspaces are created when the pixel calibration is read from the calibration database.
        # barscan_GPSANS_detector1_20200103:	0.393216 MB
        # tubewidth_GPSANS_detector1_20200130:	0.393216 MB
        DeleteWorkspace("barscan_GPSANS_detector1_20200103")
        DeleteWorkspace("tubewidth_GPSANS_detector1_20200130")


if __name__ == "__main__":
    pytest.main([__file__])
