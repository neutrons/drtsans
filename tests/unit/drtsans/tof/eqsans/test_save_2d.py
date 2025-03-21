import tempfile
import pytest
from os.path import join
from drtsans.save_2d import save_nist_dat, save_nexus
from mantid.simpleapi import LoadNexus, CompareWorkspaces
import numpy as np


@pytest.mark.datarepo
def test_save_nist_dat(datarepo_dir, clean_workspace):
    filename = join(datarepo_dir.eqsans, "test_save_output/EQSANS_68200_Iqxy.nxs")
    reference_filename = join(datarepo_dir.eqsans, "test_save_output/EQSANS_68200_Iqxy.dat")
    ws = LoadNexus(filename)
    clean_workspace(ws)
    with tempfile.NamedTemporaryFile("r+") as tmp:
        save_nist_dat(ws, tmp.name)
        output = np.loadtxt(
            tmp.name,
            dtype={"names": ("Qx", "Qy", "I", "dI"), "formats": ("f", "f", "f", "f")},
            skiprows=2,
        )
        reference = np.loadtxt(
            reference_filename,
            dtype={"names": ("Qx", "Qy", "I", "dI"), "formats": ("f", "f", "f", "f")},
            skiprows=2,
        )
        assert np.allclose(output["Qx"], reference["Qx"], atol=1e-6)
        assert np.allclose(output["Qy"], reference["Qy"], atol=1e-6)
        assert np.allclose(output["I"], reference["I"], atol=1e-6)
        assert np.allclose(output["dI"], reference["dI"], atol=1e-6)


@pytest.mark.datarepo
def test_save_nexus(datarepo_dir, clean_workspace):
    filename = join(datarepo_dir.eqsans, "test_save_output/EQSANS_68200_Iqxy.nxs")
    ws = LoadNexus(filename)
    clean_workspace(ws)
    with tempfile.NamedTemporaryFile("r+") as tmp:
        save_nexus(ws, "EQSANS 68200", tmp.name)
        output_ws = LoadNexus(tmp.name)
        clean_workspace(output_ws)
        result, messages = CompareWorkspaces(Workspace1=ws, Workspace2=output_ws)
        clean_workspace(messages)
        assert result
