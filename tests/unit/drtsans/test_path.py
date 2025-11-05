import os
import pytest
import pathlib
import stat
import sys
from tempfile import gettempdir, NamedTemporaryFile
from unittest.mock import patch

from mantid.simpleapi import CreateWorkspace
from mantid.kernel import amend_config

from drtsans.path import abspath, add_to_sys_path, exists, registered_workspace, allow_overwrite

SEARCH_ON = {}
SEARCH_OFF = {"datasearch.searcharchive": "off"}

IPTS_23732 = "/SNS/EQSANS/IPTS-23732/nexus/"
IPTS_19800 = "/SNS/EQSANS/IPTS-19800/nexus/"
IPTS_22699 = "/HFIR/CG3/IPTS-22699/nexus/"


@pytest.mark.mount_eqsans
@pytest.mark.parametrize(
    "hint, fullpath",
    [
        ("EQSANS_106026", IPTS_23732 + "EQSANS_106026.nxs.h5"),
        ("EQSANS106027", IPTS_23732 + "EQSANS_106027.nxs.h5"),
        ("EQSANS_88974.nxs.h5", IPTS_19800 + "EQSANS_88974.nxs.h5"),
    ],
    ids=("EQSANS_106026", "EQSANS_106026", "EQSANS_88974"),
)
def test_abspath_with_archivesearch(hint, fullpath, has_sns_mount):
    """Test files that require archive search in ONCat

    Note: ONCat returns the path on the analysis cluster, and Mantid's FileFinder checks that the path exists before
    returning it. This test therefore requires the SNS mount.
    """
    if not has_sns_mount:
        pytest.skip("Do not have /SNS properly mounted on this system")

    assert abspath(hint, search_archive=True) == fullpath


@pytest.mark.parametrize(
    "hint",
    ["randomname", "EQSANS_906026"],
    ids=("randomname", "EQSANS_906026"),
)
def test_abspath_file_not_exist(hint):
    with pytest.raises(RuntimeError):
        found = abspath(hint, search_archive=False)
        assert False, 'found "{}" at "{}"'.format(hint, found)


@pytest.mark.parametrize(
    "hint, instr, ipts, fullpath",
    [
        ("EQSANS_106026", "", 23732, IPTS_23732 + "EQSANS_106026.nxs.h5"),
        ("EQSANS106027", "", 23732, IPTS_23732 + "EQSANS_106027.nxs.h5"),
        ("EQSANS_88974.nxs.h5", "", 19800, IPTS_19800 + "EQSANS_88974.nxs.h5"),
        ("5709", "CG3", 22699, IPTS_22699 + "CG3_5709.nxs.h5"),
        ("5709", "CG3", 24740, IPTS_22699 + "CG3_5709.nxs.h5"),  # wrong proposal
    ],
    ids=(
        "EQSANS_106026",
        "EQSANS_106026",
        "EQSANS_88974",
        "CG3_5709",
        "CG3_5709_bad_proposal",
    ),
)
def test_abspath_with_ipts(hint, instr, ipts, fullpath):
    """Test constructing the path from file hint, instrument and IPTS"""

    def is_fullpath(path):
        return path == fullpath

    with (
        amend_config(data_dir=[], data_dir_insert_mode="replace"),
        patch("os.path.exists") as mock_path_exists,
        patch("drtsans.path.FileFinder") as mock_file_finder,
    ):
        mock_path_exists.side_effect = is_fullpath
        mock_file_finder.getFullPath.return_value = fullpath
        assert abspath(hint, instrument=instr, ipts=ipts) == fullpath


@pytest.mark.datarepo
def test_abspath_with_directory(datarepo_dir):
    filename = os.path.join(datarepo_dir.biosans, "CG3_5705.nxs.h5")
    assert abspath("CG3_5705", directory=datarepo_dir.biosans, search_archive=False) == filename
    assert (
        abspath(
            "5705",
            instrument="CG3",
            directory=datarepo_dir.biosans,
            search_archive=False,
        )
        == filename
    )


@pytest.mark.mount_eqsans
@pytest.mark.parametrize(
    "hint, found",
    [("EQSANS_106026", True), ("EQSANS106027", True), ("EQSANS_88974.nxs.h5", True)],
)
def test_exists_with_archivesearch(hint, found, reference_dir, has_sns_mount):
    """Test files that require archive search in ONCat

    Note: ONCat returns the path on the analysis cluster, and Mantid's FileFinder checks that the path exists before
    returning it. This test therefore requires the SNS mount.
    """
    if not has_sns_mount:
        pytest.skip("Do not have /SNS properly mounted on this system")

    with amend_config(SEARCH_ON, data_dir=reference_dir.eqsans):
        assert exists(hint) == found  # allows verifying against True and False


@pytest.mark.datarepo
@pytest.mark.parametrize(
    "hint, found",
    [("EQSANS_105428", True), ("EQSANS105428", True), ("EQSANS_88975.nxs.h5", True)],
)
def test_exists_without_archivesearch(hint, found, datarepo_dir):
    with amend_config(data_dir=datarepo_dir.eqsans, **SEARCH_OFF):
        assert exists(hint) == found  # allows verifying against True and False


def test_registered_workspace(temp_workspace_name):
    w_name = temp_workspace_name()
    assert registered_workspace(w_name) is False
    w = CreateWorkspace(DataX=[1], Datay=[1], OutputWorkspace=w_name)
    assert registered_workspace(w_name) is True
    assert registered_workspace(w) is True


def test_allow_overwrite(cleanfile):
    tmpdir = gettempdir()
    # create an empty file
    tmpfile = NamedTemporaryFile(dir=tmpdir, delete=False)
    tmpfile.close()
    cleanfile(tmpfile.name)  # remove the file when test finishes

    # check if others write permission is false
    path = pathlib.Path(tmpfile.name)
    assert not bool(path.stat().st_mode & stat.S_IWOTH)
    allow_overwrite(tmpdir)
    # check permissions
    assert bool(path.stat().st_mode & stat.S_IWUSR), "user writable"
    assert bool(path.stat().st_mode & stat.S_IWGRP), "group writable"
    assert bool(path.stat().st_mode & stat.S_IWOTH), "world writable"
    # delete file


class TestAddToSysPath:
    """Test suite for the add_to_sys_path context manager."""

    def test_path_added_to_sys_path(self, tmp_path):
        """Test that the path is added to sys.path."""
        test_path = str(tmp_path)
        assert test_path not in sys.path
        with add_to_sys_path(test_path):
            assert test_path in sys.path
            assert sys.path[0] == test_path  # Should be first
        assert test_path not in sys.path  # path removed after context

    def test_path_removed_on_exception(self, tmp_path):
        """Test that the path is removed even if an exception occurs."""
        test_path = str(tmp_path)
        with pytest.raises(ValueError):
            with add_to_sys_path(test_path):
                assert test_path in sys.path
                raise ValueError("Test exception")
        assert test_path not in sys.path

    def test_reduce_EQSANS_module_removed_by_default(self, tmp_path):
        """Test that reduce_EQSANS module is removed from sys.modules by default."""
        test_path = str(tmp_path)
        sys.modules["reduce_EQSANS"] = object()  # Mock module
        with add_to_sys_path(test_path):  # default clean_module=True
            assert "reduce_EQSANS" not in sys.modules

    def test_reduce_EQSANS_module_preserved_when_disabled(self, tmp_path):
        """Test that reduce_EQSANS module is preserved when clean_module=False."""
        test_path = str(tmp_path)
        mock_module = object()
        sys.modules["reduce_EQSANS"] = mock_module

        with add_to_sys_path(test_path, clean_module=False):
            assert sys.modules["reduce_EQSANS"] is mock_module
        del sys.modules["reduce_EQSANS"]  # Cleanup

    def test_reduce_EQSANS_not_present_no_error(self, tmp_path):
        """Test that no error occurs if reduce_EQSANS is not in sys.modules."""
        test_path = str(tmp_path)
        if "reduce_EQSANS" in sys.modules:
            del sys.modules["reduce_EQSANS"]
        # Should not raise any exception
        with add_to_sys_path(test_path):
            pass

    def test_original_sys_path_preserved(self, tmp_path):
        """Test that the original sys.path is preserved after context exit."""
        test_path = str(tmp_path)
        original_path = sys.path.copy()
        with add_to_sys_path(test_path):
            pass
        assert sys.path == original_path

    def test_path_added_at_beginning(self, tmp_path):
        """Test that the path is inserted at index 0 (highest priority)."""
        test_path = str(tmp_path)
        original_first = sys.path[0] if sys.path else None
        with add_to_sys_path(test_path):
            assert sys.path[0] == test_path
            if original_first:
                assert sys.path[1] == original_first


if __name__ == "__main__":
    pytest.main([__file__])
