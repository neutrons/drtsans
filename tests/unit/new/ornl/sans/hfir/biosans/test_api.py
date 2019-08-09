from mantid import mtd
from ornl.sans.hfir.biosans import load


def test_api_load(biosans_f):
    ws = load(filename=biosans_f['beamcenter'])
    assert ws.name() == "BioSANS_exp402_scan0006_0001"

    ws_name = "xpto"
    ws = load(filename=biosans_f['beamcenter'], output_workspace=ws_name)
    assert ws.name() == ws_name
    assert ws_name in mtd.getObjectNames()
