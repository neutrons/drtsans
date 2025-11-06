"""
Script digested by Mantid algorithm LoadLiveData.
LoadLiveData creates child algorithm "RunPythonScript" and runs it in a separate python interpreter process as
    RunPythonScript(InputWorkspace=input, Filename="reduce_EQSANS_live_post_proc.py")
where "input" is the EventWorkspace containing the events accumulated up to the time when the script is run.
"""

from contextlib import contextmanager
import io
import logging
import os
from os import makedirs
from shutil import copytree
import tempfile

from mantid.simpleapi import LoadEmptyInstrument
from mantid.dataobjects import EventWorkspace

from drtsans.path import add_to_sys_path
from drtsans.samplelogs import SampleLogs

LOG_NAME = "livereduce"  # same as the python logger used by service "livereduce"
LOG_FILE = "/var/log/SNS_applications/livereduce.log"
GLOBAL_AR_DIR = "/SNS/EQSANS/shared/autoreduce"


def events_file_exists(events: EventWorkspace) -> bool:
    ipts = SampleLogs(events).experiment_identifier.value[5:]  # e.g. "12345" when having IPTS-12345
    run_number = str(events.getRunNumber())  # e.g. "105584"
    events_file_path = f"/SNS/EQSANS/IPTS-{ipts}/nexus/EQSANS_{run_number}.nxs.h5"
    return os.path.isfile(events_file_path)


@contextmanager
def configure_error_buffer():
    """Context manager for error log buffer and handler."""
    error_log_buffer = io.StringIO()
    error_log_handler = logging.StreamHandler(error_log_buffer)
    error_log_handler.setLevel(logging.ERROR)
    error_log_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logging.getLogger().addHandler(error_log_handler)

    try:
        yield error_log_buffer
    finally:
        logging.getLogger().removeHandler(error_log_handler)
        error_log_handler.close()
        error_log_buffer.close()


def livereduce(events: EventWorkspace, publish=True):
    # instantiate logging components
    logger = logging.getLogger(LOG_NAME)
    logger.info("Starting EQSANS live reduction post-processing")

    # final output directory
    run = str(events.getRunNumber())
    ipts = SampleLogs(events).experiment_identifier.value[5:]  # e.g. "12345" when having IPTS-12345
    output_dir = f"/SNS/EQSANS/IPTS-{ipts}/shared/autoreduce/{run}"

    with add_to_sys_path(GLOBAL_AR_DIR):  # "/SNS/EQSANS/shared/autoreduce"
        # imports from the autoreduction script reduce_EQSANS.py
        from reduce_EQSANS import footer, LogContext, reduce_events, save_report, upload_report

        with configure_error_buffer() as error_buffer:
            log_context = LogContext(logger=logger, logfile=LOG_FILE, error_buffer=error_buffer)
            with tempfile.TemporaryDirectory(prefix="livereduce_") as temp_dir:
                report = reduce_events(events, temp_dir, log_context)
                report += footer(events, output_dir, log_context)  # notice we pass output_dir here
                save_report(report, os.path.join(temp_dir, f"EQSANS_{run}.html"), logger)  # save to disk
                if events_file_exists(events):
                    logger.warning("Translation service produced the event file for the run during live reduction.")
                else:
                    if publish:
                        upload_report(run, report, logger)  # upload to live data server
                    # copy all output files from temp_dir to output_dir
                    makedirs(output_dir, exist_ok=True)  # create a parent directory if it does not exist
                    copytree(temp_dir, output_dir, dirs_exist_ok=True)


if __name__ == "__main__":
    if events_file_exists(input):
        logging.getLogger(LOG_NAME).info("Translation service produced event file for the run. Do not live-reduce.")
    else:
        livereduce(input, publish=True)
    # Algorithm StartLiveData as invoked by livereduce.py requires an output workspace of name "result"
    output = LoadEmptyInstrument(InstrumentName="EQSANS", OutputWorkspace="result")
