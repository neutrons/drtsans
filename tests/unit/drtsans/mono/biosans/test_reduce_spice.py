from unittest import mock

import pytest

from drtsans.mono.biosans.reduce_spice import reduce_biosans_nexus


class _StopReduction(Exception):
    """Raised by a mock to stop reduce_biosans_nexus once its configuration is built."""


@pytest.mark.parametrize(
    "keywords, expected",
    [
        ({}, (0.2, 0.5)),
        ({"auto_wedge_contrast_min": 0.3, "auto_wedge_center_tolerance": 0.8}, (0.3, 0.8)),
    ],
)
def test_spice_reduction_auto_wedge_options(keywords, expected, tmp_path):
    """Test that reduce_biosans_nexus writes the auto-wedge acceptance options into the configuration."""
    module = "drtsans.mono.biosans.reduce_spice"
    pair, q_range = (1, 1), [0.003, 0.04]
    arguments = [17240, 318, ["sample"], [pair], [None], [pair], [None], pair, pair, ["0.1"], "", True, pair, pair]
    arguments += [str(tmp_path), "main.nxs", "wing.nxs", "4.05e-9", None, "", "", 100, 100, "wedge", "log", 25, 25]
    arguments += [q_range, q_range, q_range, False, None, None, 0.003, 0.04, 0.01, 0.5, 1.0, 1.0, 2.0]
    arguments += [q_range] * 6 + [25, None]
    with (
        mock.patch(f"{module}.map_to_nexus", side_effect=lambda *args, **_kwargs: ["CG3_1.nxs.h5"] * len(args[3])),
        mock.patch(f"{module}.reduction_parameters", side_effect=_StopReduction) as mock_parameters,
    ):
        with pytest.raises(_StopReduction):
            reduce_biosans_nexus(*arguments, **keywords)

    configuration = mock_parameters.call_args.args[0]["configuration"]
    assert configuration["autoWedgeSignalToNoiseMin"] == 2.0
    assert (configuration["autoWedgeContrastMin"], configuration["autoWedgeCenterTolerance"]) == expected
