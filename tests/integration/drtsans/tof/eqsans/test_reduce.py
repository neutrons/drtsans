import pytest
from os.path import join as pj
from mantid.simpleapi import LoadNexus, SumSpectra, CompareWorkspaces

from drtsans.tof.eqsans import reduce
from drtsans.settings import amend_config, unique_workspace_dundername as uwd


@pytest.mark.datarepo
def test_load_w(datarepo_dir, clean_workspace, temp_workspace_name):
    with amend_config(data_dir=datarepo_dir.eqsans):
        _w0 = reduce.load_w(
            "EQSANS_92353",
            output_workspace=uwd(),
            low_tof_clip=500,
            high_tof_clip=2000,
            dw=0.1,
        )
        _w1 = SumSpectra(_w0, OutputWorkspace=_w0.name())
        clean_workspace(_w1)
        fn = pj(datarepo_dir.eqsans, "test_reduce", "compare", "ref_load_w.nxs")
        _w2 = LoadNexus(fn, OutputWorkspace=temp_workspace_name())
        assert CompareWorkspaces(_w1, _w2)


if __name__ == "__main__":
    pytest.main([__file__])
