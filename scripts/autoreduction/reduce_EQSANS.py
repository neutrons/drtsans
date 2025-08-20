#!/usr/bin/env python
"""Autoreduction script for EQSANS"""

import argparse
import logging
import requests
import os
import warnings

from finddata.publish_plot import plot_heatmap, publish_plot
from mantid.simpleapi import Integration, LoadEventNexus, mtd

import numpy as np

# silently ignore all types of numerical errors (like divide by zero, overflow, etc.)
np.seterr(all="ignore")
warnings.filterwarnings("ignore", module="numpy")
CONDA_ENV = "sans"

TUBES_PER_EIGHTPACK = 8
TUBES_IN_DETECTOR1 = 192
PIXELS_PER_TUBE = 256


def reduce_events_file(events_file: str) -> tuple:
    """
    Process a Nexus events file to extract and transform detector data.

    Parameters
    ----------
    events_file
        Path to the Nexus events file to be processed.

    Returns
    -------
    tuple
        A tuple containing:
        - run_number (int): The run number extracted from the events file.
        - x (numpy.ndarray): Array of tube indices.
        - y (numpy.ndarray): Array of pixel indices.
        - z (numpy.ma.MaskedArray): Log-transformed and masked intensity data.

    Notes
    -----
    - If the specified file does not exist, the function will terminate the program.
    - The function uses Mantid's `LoadEventNexus` and `Integration` to process the data.
    """
    if not os.path.isfile(events_file):
        raise FileNotFoundError(f"data file {events_file} not found")

    events = LoadEventNexus(Filename=events_file, outputWorkspace=mtd.unique_hidden_name())
    intensities = Integration(InputWorkspace=events, outputWorkspace=mtd.unique_hidden_name())
    data = intensities.extractY().reshape(-1, TUBES_PER_EIGHTPACK, PIXELS_PER_TUBE).T
    data2 = data[:, [0, 4, 1, 5, 2, 6, 3, 7], :]  # tube indexes within an eightpack
    data2 = data2.transpose().reshape(-1, PIXELS_PER_TUBE)
    z = np.ma.masked_where(data2 < 1, data2)

    x = np.arange(TUBES_IN_DETECTOR1) + 1
    y = np.arange(PIXELS_PER_TUBE) + 1
    z = np.log(np.transpose(z))
    return events.getRunNumber(), x, y, z


def upload_plot(run_number: str, plot_div: str):
    """
    Upload a plot to the livedata server for the EQSANS instrument.

    Parameters
    ----------
    run_number
        The run number associated with the plot.
    plot_div
        The HTML div containing the plot to be uploaded.

    Notes
    -----
    - The function uses the `publish_plot` method to upload the plot.
    - If the upload fails due to an HTTP error, the exception is logged.
    """
    try:
        publish_plot("EQSANS", run_number, files={"file": plot_div})
    except requests.HTTPError as e:
        logging.exception(f"Publish plot failed with error {e}")


def save_plot(plot_div: str, report_file: str):
    """
    Save an HTML plot to a file with the necessary JavaScript engine for rendering.

    Parameters
    ----------
    plot_div
        The HTML div containing the plot to be saved.
    report_file
        The path to the file where the HTML report will be saved.

    Notes
    -----
    - The saved file includes the Plotly JavaScript library to ensure the plot can be displayed in a web browser.
    - The function wraps the provided `plot_div` with the required HTML structure.
    """
    # add the javascript engine so that the report can be displayed in a web browser
    prefix = """<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Plotly Chart</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>

    """
    suffix = """

    </body>
    </html>
    """
    with open(report_file, "w") as f:
        f.write(prefix + plot_div + suffix)


def parse_command_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the EQSANS autoreduction script.

    Returns
    -------
    argparse.Namespace
        An object containing the parsed command-line arguments:
        - events_file (str): Path to the Nexus events file.
        - outdir (str): Output directory path.
        - report_file (str, optional): Path to save the HTML report file.
        - no_publish (bool): Flag to disable uploading the HTML report to the livedata server.

    Notes
    -----
    - The `--report_file` argument allows specifying a file name or a full path.
      If only a file name is provided, the file is saved in the output directory.
    - The `--no_publish` flag prevents the report from being uploaded to the server.
    """
    parser = argparse.ArgumentParser(description="Autoreduction script for EQSANS")
    parser.add_argument("events_file", type=str, help="Path to the Nexus events file.")
    parser.add_argument("outdir", type=str, help="Output directory path.")
    parser.add_argument(
        "--report_file",
        type=str,
        help="Save the HTML report to file. If only the file name is given,the file is saved in the output directory",
    )
    parser.add_argument("--no_publish", action="store_true", help="Do not upload HTML report to the livedata server.")
    return parser.parse_args()


def main():
    args = parse_command_arguments()

    # Load and process the events file to extract data
    run_number, x, y, z = reduce_events_file(args.events_file)

    # Generate the plot as an HTML div
    plot_div = plot_heatmap(
        run_number,
        x.tolist(),
        y.tolist(),
        z.tolist(),
        x_title="Tube",
        y_title="Pixel",
        x_log=False,
        y_log=False,
        instrument="EQSANS",
        publish=False,
    )

    # Upload the plot to the livedata server if requested
    if not args.no_publish:
        upload_plot(run_number, plot_div)

    # Save the plot to a file if requested
    if args.report_file:
        # Find if we should save to the output directory
        file_path = args.report_file
        if file_path and os.path.basename(file_path) == file_path:
            file_path = os.path.join(args.outdir, file_path)
        save_plot(plot_div, file_path)


if __name__ == "__main__":
    main()
