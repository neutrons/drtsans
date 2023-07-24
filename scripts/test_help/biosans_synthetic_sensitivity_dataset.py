"""
The purpose of this script is to load a set of Nexus Events files and:
- integrate all events per pixel, yielding a Workspace2D with one bin per pixel and units of time-of-flight
- insert the monitor counts in log "monitor"
- clone the intensities of the wing detector onto the midrange detector

The Nexus Event files corresponding to runs 4835, 4830, and 4831 are used in the integration test for the
sensitivity of the Wing detector. The Files produced by this script can be used to test the sensitivity of
the Midrange Detector. See Fixture ``biosans_synthetic_sensitivity_dataset`` for detailed use

"""
# local imports
from drtsans.load import __monitor_counts
from drtsans.mono.biosans.simulated_intensities import clone_component_intensities, insert_midrange_detector
from drtsans.samplelogs import SampleLogs
from drtsans.settings import unique_workspace_dundername

# third-party imports
from mantid.simpleapi import Integration, LoadEventNexus, SaveNexus

# standard library imports
import os


def clone_intensities_to_midrange(filepath_in, filepath_out):
    monitor_count = __monitor_counts(filepath_in)
    input_workspace = unique_workspace_dundername()
    LoadEventNexus(Filename=filepath_in, OutputWorkspace=input_workspace)
    Integration(InputWorkspace=input_workspace, OutputWorkspace=input_workspace)
    SampleLogs(input_workspace).insert("monitor", monitor_count)
    output_workspace = unique_workspace_dundername()
    insert_midrange_detector(input_workspace, output_workspace=output_workspace)
    clone_component_intensities(
        output_workspace, input_component="wing_detector", output_component="midrange_detector"
    )
    SaveNexus(output_workspace, filepath_out)


if __name__ == "__main__":
    input_data_dir = "/HFIR/CG3/IPTS-23782/nexus"
    assert os.path.exists(input_data_dir), "Unable to generate the dataset."

    runs = dict(FLOOD_RUNS=4835, DIRECT_BEAM_RUNS=4830, TRANSMISSION_RUNS=4831, TRANSMISSION_FLOOD_RUNS=4835)
    # output_data_dir = "/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/biosans/synthetic_sensitivity"
    output_data_dir = "/tmp"
    for run in runs.values():
        clone_intensities_to_midrange(f"{input_data_dir}/CG3_{run}.nxs.h5", f"{output_data_dir}/CG3_{run}.nxs")
