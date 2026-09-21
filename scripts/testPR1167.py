"""Manual smoke test for PR 1167 (requires the SNS data mount and Mantid)."""

from drtsans.beam_finder import fbc_options_json
from drtsans.redparams import reduction_parameters
from drtsans.tof.eqsans import load_all_files, reduce_single_configuration


INPUT = {
    "instrumentName": "EQSANS",
    "iptsNumber": "37596",
    "outputFileName": "PR1167_188546",
    "sample": {"runNumber": "188546"},
    "beamCenter": {
        "method": "center_of_mass",
        "runNumber": "188546",
        "useFallbackBeamCenter": True,
        "fallbackBeamCenter": [0.025239, 0.0170801],
        "com_centering_options": {"IntegrationRadius": 0.1},
    },
    "configuration": {
        "outputDir": "/tmp/PR1167_188546",
        "QbinType": "linear",
        "numQBins": 60,
    },
}


config = reduction_parameters(INPUT, permissible=True)
print("Forwarded beam-center options:", fbc_options_json(config))
loaded = load_all_files(config, prefix="PR1167")
reduce_single_configuration(loaded, config, prefix="PR1167")
