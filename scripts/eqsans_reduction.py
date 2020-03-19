import json
import os
import sys
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402
from drtsans.tof.eqsans.api import load_all_files, reduce_single_configuration, plot_reduction_output  # noqa E402

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError(".!.reduction code requires a parameter json string or file")
    if os.path.isfile(sys.argv[1]):
        json_filename = sys.argv[1]
        msapi.logger.warning('...json file is used.')
        with open(sys.argv[1], 'r') as fd:
            reduction_input = json.load(fd)
    else:
        msapi.logger.warning('...json string is used.')
        json_filename = " ".join(sys.argv[1:])
        reduction_input = json.loads(json_filename)

    # chekcing if output directory exists, if it doesn't, creates the folder
    output_dir = reduction_input["configuration"]["outputDir"]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    loaded = load_all_files(reduction_input)
    msapi.logger.warning('...loading completed.')
    out = reduce_single_configuration(loaded, reduction_input)
    msapi.logger.warning('...single reduction completed.')
    plot_reduction_output(out, reduction_input)
    msapi.logger.warning('...plotting completed.')
