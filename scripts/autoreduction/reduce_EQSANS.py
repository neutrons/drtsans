#!/usr/bin/env python
"""Autoreduction script for EQSANS"""

import argparse
from copy import deepcopy
import datetime
import io
import json
import logging
import os
import requests
import time
from typing import Union
import warnings

from plot_publisher import plot_heatmap, publish_plot
from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import DeleteWorkspace, Integration, LoadEventNexus, mtd
from mantid.utils.logging import log_to_python as mtd_log_to_python
import numpy as np
import plotly.offline as pyo

import drtsans
from drtsans.redparams import reduction_parameters, update_reduction_parameters
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans.api import (
    load_all_files,
    reduce_single_configuration,
    plot_reduction_output,
    plotly_reduction_output,
)
from drtsans.tof.eqsans.meta_data import is_sample_run

# silently ignore all types of numerical errors (like divide by zero, overflow, etc.)
np.seterr(all="ignore")
warnings.filterwarnings("ignore", module="numpy")
CONDA_ENV = "sans"

LOG_NAME = "autoreduce"
AUTOREDUCE_DIR = "/SNS/EQSANS/shared/autoreduce"
AUTOREDUCE_IPTS_DIR = "/SNS/EQSANS/IPTS-{ipts}/shared/autoreduce"  # e.g. /SNS/EQSANS/IPTS-12345/shared/autoreduce

TUBES_PER_EIGHTPACK = 8
TUBES_IN_DETECTOR1 = 192
PIXELS_PER_TUBE = 256


def filelink(filepath: str) -> str:
    """Create an HTML link to a file path.

    Parameters
    ----------
    filepath : str
        The file path to be linked.

    Returns
    -------
    str
        An HTML anchor tag linking to the specified file path.
    """
    return f"<a href='file://{filepath}' target='_blank'>{filepath}</a>"


def configure_logger(output_dir: str):
    """Configure logging for the autoreduction process.

    Sets up file and error logging handlers for the autoreduction workflow. Redirects Mantid
    logging to Python's logging system and creates a buffer to capture error messages separately.

    Parameters
    ----------
    output_dir : str
        Directory path where the log file will be saved. The log file will be named 'autoreduce.log'.

    Returns
    -------
    io.StringIO
        A StringIO buffer containing error-level log messages. This buffer can be used to
        display errors in the final HTML report.

    Notes
    -----
    - Mantid logging is redirected to Python logging at INFO level
    - A file handler writes all INFO and above messages to `autoreduce.log`
    - A separate StringIO handler captures only ERROR and above messages
    - The root logger for 'autoreduce' is set to INFO level
    """
    # redirect Mantid logging to python logging
    mtd_log_to_python("information")
    logging.getLogger("Mantid").setLevel(logging.INFO)

    # create a file handler
    fileHandler = logging.FileHandler(os.path.join(output_dir, f"{LOG_NAME}.log"))
    fileHandler.setLevel(logging.INFO)
    logformat = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    fileHandler.setFormatter(logging.Formatter(logformat))
    logging.getLogger().addHandler(fileHandler)

    # Create a StringIO buffer and handler for error messages
    error_log_buffer = io.StringIO()
    error_log_handler = logging.StreamHandler(error_log_buffer)
    error_log_handler.setLevel(logging.ERROR)
    error_log_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logging.getLogger().addHandler(error_log_handler)

    # add the handlers to the python root logger

    logging.getLogger(LOG_NAME).setLevel(logging.INFO)
    return error_log_buffer


def intensity_array(events: EventWorkspace) -> tuple:
    """
    Integrate intensities of each detector pixel.

    Parameters
    ----------
    events
        Raw events to be processed.

    Returns
    -------
    tuple
        A tuple containing:
        - x (numpy.ndarray): Array of tube indices.
        - y (numpy.ndarray): Array of pixel indices.
        - z (numpy.ma.MaskedArray): Log-transformed and masked intensity data.

    Notes
    -----
    - The function uses Mantid's `Integration` to process the data.
    """
    workspace_name = mtd.unique_hidden_name()
    intensities = Integration(InputWorkspace=events, OutputWorkspace=workspace_name)
    data = intensities.extractY().reshape(-1, TUBES_PER_EIGHTPACK, PIXELS_PER_TUBE).T
    data2 = data[:, [0, 4, 1, 5, 2, 6, 3, 7], :]  # tube indexes within an eightpack
    data2 = data2.transpose().reshape(-1, PIXELS_PER_TUBE)
    z = np.ma.masked_where(data2 < 1, data2)
    x = np.arange(TUBES_IN_DETECTOR1) + 1
    y = np.arange(PIXELS_PER_TUBE) + 1
    z = np.log(np.transpose(z))
    DeleteWorkspace(workspace_name)
    return x, y, z


def upload_report(run_number: str, plot_div: str):
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
        logging.getLogger(LOG_NAME).exception(f"Publish plot failed with error {e}")


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
        logging.getLogger(LOG_NAME).error(f"Plotly.js version {js_version} not found, using version 3.0.0 instead")
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


def save_report(report: str, report_file: str):
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


def reduce_non_sample(events: EventWorkspace):
    """
    Reduce events from a non-sample run by generating a heatmap plot of pixel intensities

    Parameters
    ----------
    events : EventWorkspace
        The Mantid EventWorkspace containing the raw events to be processed.

    Returns
    -------
    report (str): An HTML div containing the generated heatmap plot of pixel intensities
    """
    # Aggregate the events within each pixel as an intensity value
    x, y, z = intensity_array(events)

    # Generate the plot as an HTML div
    report = plot_heatmap(
        events.getRunNumber(),
        x.tolist(),
        y.tolist(),
        z.tolist(),
        x_title="Tube",
        y_title="Pixel",
        x_log=False,
        y_log=False,
        instrument="EQSANS",
        title="Counts per Pixel",
        publish=False,
    )
    return report


def reduce_sample(events: EventWorkspace, output_dir: str):
    """Reduce events from a sample run and generate comprehensive output files and plots.

    This function performs the complete reduction workflow for sample runs, including loading
    reduction options, validating configuration, reducing data, and generating output files
    and plots.

    Parameters
    ----------
    events : EventWorkspace
        The Mantid EventWorkspace containing the raw sample events to be reduced.
    output_dir : str
        Base output directory path for reduction files. If this matches the IPTS-specific
        autoreduce directory, a subdirectory with the run number will be created for outputs.

    Returns
    -------
    str
        An HTML string containing the complete reduction report, including:
        - Pixel intensity heatmap (from reduce_non_sample)
        - Plotly reduction output plots (pixel integrated intensity heatmap, I(Qx,Qy) heatmap, I(Q) profiles)

    Notes
    -----
    - Searches for reduction options in order of priority:
      1. `reduction_options_{run_number}.json` in output_dir
      2. `reduction_options.json` in IPTS autoreduce directory
      3. `reduction_options.json` in shared autoreduce directory
    - Automatically amends reduction parameters with run-specific information
    - Saves reduced data files, HDF5 log, and plots (*.png) to the output directory
    - Saves the final comprehensive reduction options to `reduction_options_{run_number}.json`
    """
    # The autoreduction service will pass AUTOREDUCE_IPTS_DIR as output directory.
    # We add run number as subdir, otherwise there will be too many files under AUTOREDUCE_IPTS_DIR
    run_number = str(events.getRunNumber())  # e.g. "105584"
    ipts = SampleLogs(events).experiment_identifier.value[5:]  # e.g. "12345" when having IPTS-12345
    reduced_files_dir = output_dir
    if output_dir == AUTOREDUCE_IPTS_DIR.format(ipts=ipts):
        reduced_files_dir = os.path.join(output_dir, run_number)

    # find most appropriate reduction options and amend if necessary
    amendment = {}
    #  Example: /SNS/EQSANS/IPTS-12345/shared/autoreduce/reduction_options_105584.json
    reduction_options_path = os.path.join(output_dir, f"reduction_options_{run_number}.json")
    if os.path.exists(reduction_options_path) is False:
        #  Example: /SNS/EQSANS/IPTS-12345/shared/autoreduce/reduction_options.json
        amendment = {
            "iptsNumber": ipts,
            "sample": {"runNumber": run_number},
            "outputFileName": f"EQSANS_{run_number}",  # prefix for all output files
            "configuration": {"outputDir": reduced_files_dir},
        }
        reduction_options_path = os.path.join(AUTOREDUCE_IPTS_DIR.format(ipts=ipts), "reduction_options.json")
        if os.path.exists(reduction_options_path) is False:
            #  Fallback: /SNS/EQSANS/shared/autoreduce/reduction_options.json
            reduction_options_path = os.path.join(AUTOREDUCE_DIR, "reduction_options.json")
            amendment["beamCenter"] = {"runNumber": run_number}
    footer = (
        "<b>Note:</b> the links below can't be open directly. Instead, copy the link address "
        "and paste it in a new browser tab.<br>\n"
    )
    footer += "<table border='0'><tr><td>\n"
    footer += f"Input reduction options loaded from</td><td>{filelink(reduction_options_path)}</td></tr>\n"
    footer += "</table>\n"

    # load options and validate
    with open(reduction_options_path, "r") as f:
        raw_options = json.load(f)
    input_config = reduction_parameters(raw_options, validate=False, permissible=True)
    input_config = update_reduction_parameters(input_config, amendment, validate=True, permissible=True)
    final_input_config = deepcopy(input_config)  # input_config will be modified during reduction

    # load files and reduce
    loaded = load_all_files(input_config)
    output = reduce_single_configuration(loaded, input_config)

    # create plot images for all intensity profiles
    plot_reduction_output(output, input_config)  # save as *.PNG files

    # collect all HTML elements into a single report
    report = reduce_non_sample(events) + "<hr>\n" + plotly_reduction_output(output, input_config) + "<hr>\n"

    # Save the input reduction options
    reduction_options_path = os.path.join(output_dir, f"reduction_options_{run_number}.json")
    with open(reduction_options_path, "w") as config_file:
        json.dump(final_input_config, config_file, indent=4, sort_keys=True)

    return report + footer


def footer(events: EventWorkspace, output_dir: str) -> str:
    """Generate an HTML footer with reduction metadata and file locations.

    Creates an HTML table containing information about the reduction process, including
    software versions, file paths, timing information, and reduction parameters.

    Parameters
    ----------
    events : EventWorkspace
        The Mantid EventWorkspace containing the raw events data. Used to extract
        run number, IPTS number, proton charge statistics, and sample logs.
    output_dir : str
        Base output directory path where reduction files are saved. If this matches
        the IPTS-specific autoreduce directory, a subdirectory with the run number
        will be used for the actual output files.

    Returns
    -------
    str
        An HTML string containing a table with reduction metadata, including:
        - drtsans version and documentation links
        - Log file location
        - Proton charge duration
        - Reduction options file paths
        - Output files directory
        - Reduction timestamp

    Notes
    -----
    - The function determines the actual reduced files directory based on the output_dir
    - If output_dir matches the IPTS autoreduce pattern, files are saved to a run-specific subdirectory
    - The footer includes horizontal rules (<hr>) for visual separation in the HTML report
    """
    run_number = str(events.getRunNumber())
    ipts = SampleLogs(events).experiment_identifier.value[5:]  # e.g. "12345" when having IPTS-12345
    reduced_files_dir = output_dir
    if output_dir == AUTOREDUCE_IPTS_DIR.format(ipts=ipts):
        reduced_files_dir = os.path.join(output_dir, run_number)
    reduction_options_path = os.path.join(output_dir, f"reduction_options_{run_number}.json")
    docs = "<a href='https://drtsans.readthedocs.io/latest/index.html' target='_blank'>drtsans</a>"
    version = drtsans.__version__
    release = f"<a href='https://github.com/neutrons/drtsans/releases/tag/v{version}' target='_blank'>{version}</a>"
    footer = "<table border='0'>\n"
    footer += f"<tr><td>Reduced with </td><td>{docs} version {release}</td></tr>\n"
    log_file = os.path.join(reduced_files_dir, "autoreduce.log")
    footer += f"<tr><td>Reduction log</td><td>{filelink(log_file)}</td></tr>\n"
    proton_charge = SampleLogs(events).proton_charge
    footer += f"<tr><td>Duration from proton charge</td><td>{proton_charge.getStatistics().duration} sec</td></tr>\n"
    footer += "<tr><td>Comprehensive input reduction options saved to</td><td>"
    footer += f"{filelink(reduction_options_path)}</td></tr>\n"
    footer += f"<tr><td>Output reduction files saved to</td><td>{filelink(reduced_files_dir)}</td></tr>\n"
    footer += f"<tr><td>Date </td><td>{datetime.datetime.now().strftime('%Y-%m-%d %I:%M %p')}</td></tr>\n"
    footer += "</table>\n<hr>\n"
    return footer


def autoreduce(args: argparse.Namespace):
    """Execute the complete autoreduction workflow for EQSANS event data.

    This is the main entry point for the autoreduction process. It orchestrates loading event data,
    performing reduction (either sample or non-sample), generating reports, handling errors, and
    optionally publishing results to the livedata server.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments containing:
        - events_file (str): Path to the Nexus events file to be reduced
        - outdir (str): Base output directory for reduction files and reports
        - no_publish (bool): If True, skip uploading the report to the livedata server

    Raises
    ------
    FileNotFoundError
        If the specified events file does not exist

    Notes
    -----
    - Configures logging to capture all reduction process information and errors
    - Automatically determines if the run is a sample or non-sample run
    - For IPTS autoreduce directories, creates a run-specific subdirectory
    - Saves an HTML report containing plots, metadata, and any error messages
    - Reports total reduction time upon completion
    - Reports total reduction time upon completion
    - If reduction fails, captures the exception in the error log
    - Optionally uploads the report to the livedata server (default: enabled)
    """
    start_time = time.time()
    error_log_buffer = configure_logger(args.outdir)  # save error messages to this buffer

    # Load events file
    if not os.path.isfile(args.events_file):
        raise FileNotFoundError(f"data file {args.events_file} not found")
    events = LoadEventNexus(Filename=args.events_file, OutputWorkspace=mtd.unique_hidden_name())
    run_number = str(events.getRunNumber())

    # Output directory for reduced files, logs, and report
    reduced_files_dir = args.outdir
    ipts = SampleLogs(events).experiment_identifier.value[5:]  # e.g. "12345" when having IPTS-12345
    if args.outdir == AUTOREDUCE_IPTS_DIR.format(ipts=ipts):  # e.g. /SNS/EQSANS/IPTS-12345/shared/autoreduce/
        reduced_files_dir = os.path.join(args.outdir, run_number)  # /SNS/EQSANS/IPTS-12345/shared/autoreduce/105584

    # reduce events
    report = ""
    try:
        report += reduce_sample(events, args.outdir) if is_sample_run(events) else reduce_non_sample(events)
        report += footer(events, args.outdir)
        minutes, seconds = divmod(time.time() - start_time, 60)
        report += f"<div>Reduction completed in: <b>{int(minutes)} min {int(seconds)} sec.</b><br></div>\n"
    except Exception:
        logging.getLogger(LOG_NAME).exception("Autoreduction failed")

    # possibly add error log messages
    error_messages = error_log_buffer.getvalue()
    if error_messages:
        report += f"<div><h3>Error Messages</h3><pre>{error_messages}</pre></div></hr>\n"

    # Save report to disk as an HTML file
    save_report(report, os.path.join(reduced_files_dir, f"EQSANS_{run_number}.html"))

    #  Upload report to the livedata server if requested
    if not args.no_publish:
        upload_report(run_number, report)


if __name__ == "__main__":
    args = parse_command_arguments()
    autoreduce(args)
