import argparse
import json
import os

from drtsans.mono.gpsans import reduction_parameters
from drtsans.mono.gpsans.api import (
    load_all_files,
    reduce_single_configuration,
    plot_reduction_output,
)


def parse_command_arguments():
    """
    Parse command-line arguments for the GPSANS reduction script.

    Returns
    -------
    argparse.Namespace
        An object containing the parsed command-line arguments:
        - config (str): Path to a JSON file or a JSON string containing reduction parameters.
        - permissible (bool): Allow permissible parameters in reduction configuration.

    Notes
    -----
    - If the config argument is a path to an existing file, it will be read as a JSON file.
    - Otherwise, it will be parsed as a JSON string (supports multiple arguments joined together).
    """
    parser = argparse.ArgumentParser(
        description="GPSANS reduction script. Accepts either a JSON file path or a JSON string."
    )
    parser.add_argument(
        "config",
        type=str,
        nargs="+",
        help="Path to a JSON file containing reduction parameters, or a JSON string.",
    )
    parser.add_argument(
        "--permissible",
        action="store_true",
        default=False,
        help="raise only a warning (not an exception) if a parameter in input 'config'"
        " is not found in the instrument's schema.",
    )
    return parser.parse_args()


def reduce_gpsans_configuration(input_config, permissible=False):
    # Serve and validate all necessary reduction parameters
    config = reduction_parameters(input_config, permissible=permissible)
    # create the output directory if non-existent
    output_dir = config["configuration"]["outputDir"]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # load all files and perform reduction
    loaded = load_all_files(config)
    out = reduce_single_configuration(loaded, config)
    # create figures in the output directory
    plot_reduction_output(out, config)


if __name__ == "__main__":
    args = parse_command_arguments()

    # Join config arguments (supports both single and multiple arguments)
    config_input = " ".join(args.config)

    if os.path.isfile(config_input):
        with open(config_input, "r") as fd:
            reduction_input = json.load(fd)
    else:
        reduction_input = json.loads(config_input)

    reduce_gpsans_configuration(reduction_input, permissible=args.permissible)
