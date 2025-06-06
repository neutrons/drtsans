from mantid.simpleapi import MaskBTP
from mantid.simpleapi import DeleteWorkspace
import numpy as np
from drtsans.mono import gpsans as sans
import drtsans
import lmfit
import pytest
import json
import os


# defining 2D Gaussian fitting functions
def Gaussian2D(x1, y1, amp, sigma_x, sigma_y, theta, x0, y0):
    a = np.cos(theta) ** 2 / (2.0 * sigma_x**2) + np.sin(theta) ** 2 / (2.0 * sigma_y**2)
    b = -np.sin(2.0 * theta) / (4.0 * sigma_x**2) + np.sin(2.0 * theta) / (4.0 * sigma_y**2)
    c = np.sin(theta) ** 2 / (2.0 * sigma_x**2) + np.cos(theta) ** 2 / (2.0 * sigma_y**2)
    amplitude = amp / (np.sqrt(2.0 * np.pi) * np.sqrt(sigma_x * sigma_y))
    val = amplitude * np.exp(-(a * (x1 - x0) ** 2 + 2.0 * b * (x1 - x0) * (y1 - y0) + c * (y1 - y0) ** 2))
    return val


def test_gaussian_fit(datarepo_dir, temp_directory):
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

    flood_file = os.path.join(datarepo_dir.gpsans, "sens_c489_bar.nxs")
    # Find beam center for main detector
    # loading beam center data
    center_filename = os.path.join(datarepo_dir.gpsans, "CG2_9188.nxs.h5")
    ws = sans.load_events(center_filename, output_workspace="ws_center", pixel_calibration=pixel_calibration_filename)
    ws = sans.transform_to_wavelength(ws)
    ws = drtsans.process_uncertainties.set_init_uncertainties(ws)
    sans.solid_angle_correction(ws)
    drtsans.apply_sensitivity_correction(ws, flood_file, min_threshold=0.5, max_threshold=1.5)
    MaskBTP(ws, Pixel="1-70,186-256")

    # fitting 2D gaussian to center data
    ws_size = ws.getNumberHistograms()
    x = np.empty(ws_size)
    y = np.empty(ws_size)
    intes = np.empty(ws_size)
    intes_err = np.empty(ws_size)
    keep = np.empty(ws_size, dtype=np.bool_)

    for i, si in enumerate(ws.spectrumInfo()):
        pos = si.position
        x[i] = pos.X()
        y[i] = pos.Y()
        keep[i] = not si.isMasked and np.isfinite(ws.readY(i)[0])
        intes[i] = ws.readY(i)[0]
        intes_err[i] = ws.readE(i)[0]

    x = x[keep]
    y = y[keep]
    intes = intes[keep]
    intes_err = intes_err[keep]

    model = lmfit.Model(
        Gaussian2D,
        independent_vars=["x1", "y1"],
        param_names=["amp", "sigma_x", "sigma_y", "theta", "x0", "y0"],
    )

    params = lmfit.Parameters()
    params.add("amp", value=ws.extractY().max())
    params.add("sigma_x", value=0.01, min=np.finfo(float).eps)  # width in x
    params.add("sigma_y", value=0.01, min=np.finfo(float).eps)  # width in y
    params.add("theta", value=np.pi / 2)
    params.add("x0", value=0.0)
    params.add("y0", value=0.0)
    params["theta"].vary = False
    results = model.fit(intes, x1=x, y1=y, weights=1.0 / intes_err, params=params)

    x0, y0, fit_results = sans.find_beam_center(
        ws,
        method="gaussian",
        centering_options={"theta": {"value": np.pi / 2.0, "vary": False}},
        solid_angle_method=None,
    )
    assert x0 == pytest.approx(results.params["x0"].value)
    assert y0 == pytest.approx(results.params["y0"].value)
    # update ref val after upgrading mantid (v6.11 -> v6.12)
    assert fit_results["amp"]["value"] == pytest.approx(720510910)
    assert fit_results["sigma_x"]["value"] == pytest.approx(0.0133174719397)
    assert fit_results["sigma_y"]["value"] == pytest.approx(0.0075641382315)
    assert fit_results["theta"]["value"] == pytest.approx(np.pi / 2.0)

    params["theta"].value = 0.0
    results = model.fit(intes, x1=x, y1=y, weights=1.0 / intes_err, params=params)
    x0, y0, fit_results = sans.find_beam_center(
        ws,
        method="gaussian",
        centering_options={"theta": {"value": 0.0, "vary": False}},
        solid_angle_method=None,
    )
    assert x0 == pytest.approx(results.params["x0"].value)
    assert y0 == pytest.approx(results.params["y0"].value)
    # update ref val after upgrading mantid (v6.11 -> v6.12)
    assert fit_results["amp"]["value"] == pytest.approx(720510910)
    assert fit_results["sigma_x"]["value"] == pytest.approx(0.0075641382315)
    assert fit_results["sigma_y"]["value"] == pytest.approx(0.0133174719397)
    assert fit_results["theta"]["value"] == pytest.approx(0.0)

    # cleanup
    DeleteWorkspace(ws)
    # NOTE:
    # the sensitivity calculation has two leftover workspaces, and it is unclear if we where
    # they were generated:
    # barscan_GPSANS_detector1_20200818:	0.393216 MB
    # tubewidth_GPSANS_detector1_20200612:	0.393216 MB
    # for testing purposes, we will manually delete them here
    DeleteWorkspace("barscan_GPSANS_detector1_20200103")
    DeleteWorkspace("tubewidth_GPSANS_detector1_20200130")


if __name__ == "__main__":
    pytest.main([__file__])
