#!/usr/bin/env python
"""Autoreduction script for EQSANS"""

import argparse
import logging
import os
import requests
from typing import Union
import warnings

from finddata.publish_plot import plot_heatmap, publish_plot
import plotly.offline as pyo
from mantid.simpleapi import Integration, LoadEventNexus, logger, mtd

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


def html_wrapper(report: Union[str, None]) -> str:
    """Wraps a report (set of <dvi> elements) in a complete HTML document

    Adds the javascript engine (PlotLy.js) address, HTML head, and body tags.

    Parameters
    ----------
    report : str
        The HTML content to be wrapped. This should contain one or more <div> elements and possibly
        summary <table> elements.

    Returns
    -------
    str
        A complete HTML document as a string, with the provided report content embedded within the body.
    """
    js_version = pyo.get_plotlyjs_version()
    url = f"https://cdn.plot.ly/plotly-{js_version}.js"
    try:
        response = requests.head(url, timeout=5)
        assert response.status_code == 200
    except (requests.RequestException, AssertionError):
        logger.error(f"Plotly.js version {js_version} not found, using version 3.0.0 instead")
        url = "https://cdn.plot.ly/plotly-3.0.0.js"

    prefix = f"""<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Plotly Chart</title>
        <script src="{url}"></script>
    </head>
    <body>

    """
    suffix = """

    </body>
    </html>
    """
    report = "" if report is None else report
    return prefix + report + suffix  # allow for report being `None`


def save_plot(report: str, report_file: str):
    """
    Save an HTML plot to a file with the necessary JavaScript engine for rendering.

    Parameters
    ----------
    report
        The HTML content to be saved. This may contain one or more <div> elements
        and possibly summary <table> elements.
    report_file
        The path to the file where the HTML report will be saved.

    Notes
    -----
    - The saved file includes the Plotly JavaScript library to ensure the plot can be displayed in a web browser.
    - The function wraps the provided `report` with the required HTML structure.
    """
    with open(report_file, "w") as f:
        f.write(html_wrapper(report))


def parse_command_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the EQSANS autoreduction script.

    Returns
    -------
    argparse.Namespace
        An object containing the parsed command-line arguments:
        - events_file (str): Path to the Nexus events file.
        - outdir (str): Output directory path.
        - no_publish (bool): Flag to disable uploading the HTML report to the livedata server.

    Notes
    -----
    - The `--no_publish` flag prevents the report from being uploaded to the server.
    """
    parser = argparse.ArgumentParser(description="Autoreduction script for EQSANS")
    parser.add_argument("events_file", type=str, help="Path to the Nexus events file.")
    parser.add_argument("outdir", type=str, help="Output directory path.")
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

    #  Save the plot to disk and upload the plot to the livedata server if requested
    os.makedirs(args.outdir, exist_ok=True)
    save_plot(plot_div, os.path.join(args.outdir, f"EQSANS_{run_number}.html"))
    if not args.no_publish:
        upload_plot(run_number, plot_div)


if __name__ == "__main__":
    main()
