import sys
import os
import json
from drtsans.mono.gpsans import load_all_files, reduce_single_configuration, plot_reduction_output

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("reduction code requires a parameter json string")
    if os.path.isfile(sys.argv[1]):
        print(sys.argv[1])
        with open(sys.argv[1], 'r') as fd:
            reduction_input = json.load(fd)
    else:
        json_string = " ".join(sys.argv[1:])
        reduction_input = json.loads(json_string)

    # chekcing if output directory exists, if it doesn't, creates the folder
    output_dir = reduction_input["configuration"]["outputDir"]
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    loaded = load_all_files(reduction_input)
    out = reduce_single_configuration(loaded, reduction_input)
    plot_reduction_output(out, reduction_input)
