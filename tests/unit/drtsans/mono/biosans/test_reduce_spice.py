from unittest import mock

from mantid.kernel import config
from mantid.simpleapi import CreateSingleValuedWorkspace, GroupWorkspaces, mtd
import pytest

from drtsans.mono.biosans.reduce_spice import all_workspace_names, clear_buffer, reduce_biosans_nexus


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


@pytest.mark.parametrize("show_hidden", ["0", "1"])
def test_clear_buffer(clean_workspace, show_hidden):
    """clear_buffer deletes visible and hidden workspaces, except those in ``keep``, whatever the Mantid setting"""
    setting = "MantidOptions.InvisibleWorkspaces"
    original = config[setting]
    config[setting] = show_hidden
    try:
        workspaces_before = set(all_workspace_names())  # may belong to other tests running in this process
        names = ["visible", "kept", "member1", "member2", "__hidden", "__kept_hidden"]
        for name in names:
            CreateSingleValuedWorkspace(1.0, OutputWorkspace=clean_workspace(name))
        GroupWorkspaces(InputWorkspaces="member1,member2", OutputWorkspace=clean_workspace("group"))
        assert {"__hidden", "__kept_hidden"} <= set(all_workspace_names())

        clear_buffer(keep=workspaces_before | {"kept", "__kept_hidden"})
        remaining = {name for name in names + ["group"] if mtd.doesExist(name)}
        assert remaining == {"kept", "__kept_hidden"}
        assert all(mtd.doesExist(name) for name in workspaces_before)
        assert config[setting] == show_hidden
    finally:
        config[setting] = original


def test_clear_buffer_without_keep():
    """Without ``keep``, clear_buffer deletes every workspace still registered, hidden ones included"""
    module = "drtsans.mono.biosans.reduce_spice"
    with mock.patch(f"{module}.mtd") as mock_mtd, mock.patch(f"{module}.all_workspace_names") as mock_names:
        mock_names.return_value = ["__hidden", "group", "member1", "visible"]
        mock_mtd.doesExist.side_effect = lambda name: name != "member1"  # deleted along with "group"
        clear_buffer()
    assert [call.args[0] for call in mock_mtd.remove.call_args_list] == ["__hidden", "group", "visible"]
