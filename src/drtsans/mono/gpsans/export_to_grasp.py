# Data exporter for importing into Grasp.

import drtsans
from drtsans.mono import gpsans as sans
from mantid.simpleapi import mtd
from drtsans.tubecollection import TubeCollection
from drtsans.geometry import panel_names
import numpy as np

import argparse
import os
import re
import time


def main():
    parser = argparse.ArgumentParser(description="drtsans GRASP")

    # Add positional arguments
    parser.add_argument("datafile", type=str, help="The datafile to be processed.")
    parser.add_argument(
        "output_directory", type=str, help="The output directory that the GRASP file will be saved to."
    )

    args = parser.parse_args()

    filename = args.datafile
    outdir = args.output_directory

    # print(filename)
    # print(outdir)
    # exit()
    # USER Input here with scan numbers etc.
    # ipts='30979'
    # run_min=74849
    # run_max=75000#36649

    # centers=['32422']

    mag_type = "magG"  # magH or magG

    output_directory = outdir

    ##########################Local contact adds  flood file location####################
    flood_file = "/HFIR/CG2/shared/drt_sensitivity/sens_c511_bar.nxs"
    # samples=[str(x) for x in list(range(run_min,run_max+1,1))]
    # sample_names=samples
    ###############Reduction start here, no user input needed below ########################

    # chekcing if output directory exists, if it doesn't, creates the folder

    output_dir = output_directory
    """
    for subfolder in ['1D','2D']:
        output_folder=os.path.join(output_dir,subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
    """

    # get_ipython().run_line_magic('matplotlib', 'inline')

    # Find beam center

    # center_filename = f"/HFIR/CG2/IPTS-{ipts}/nexus/CG2_{centers[0]}.nxs.h5"
    # ws=sans.load_events(center_filename,output_workspace='ws_center', pixel_calibration=True)
    # LoadInstrument(ws, InstrumentName='CG2', RewriteSpectraMap=False)
    # ws=sans.transform_to_wavelength(ws)
    # ws=drtsans.process_uncertainties.set_init_uncertainties(ws)
    # xc, yc, fits = sans.find_beam_center(ws)
    SAMPLE_SI_META_NAME = "CG2:CS:SampleToSi"
    SI_WINDOW_NOMINAL_DISTANCE_METER = 0

    script_start_time = time.time()
    # for idx in range(len(samples)):
    # run_filename = f"/HFIR/CG2/IPTS-{ipts}/nexus/CG2_{samples[idx]}.nxs.h5"
    run_filename = filename
    print(run_filename)
    w_f = sans.load_events_and_histogram(
        run_filename,
        sample_to_si_name=SAMPLE_SI_META_NAME,
        si_nominal_distance=SI_WINDOW_NOMINAL_DISTANCE_METER,
        output_workspace="sample_ws",
        pixel_calibration=True,
    )

    w_f = sans.transform_to_wavelength(w_f)
    w_f = drtsans.process_uncertainties.set_init_uncertainties(w_f)

    drtsans.apply_sensitivity_correction(w_f, flood_file, min_threshold=0.5, max_threshold=1.5)

    # plot_detector(w_f,backend='mpl',imshow_kwargs={'norm': LogNorm(vmin = 1)})
    run_obj = w_f.getRun()
    wavelength = run_obj["wavelength"].getStatistics().mean
    wavelength_spread = run_obj["wavelength_spread"].getStatistics().mean

    detector_dist = run_obj["sample_detector_distance"].value

    detector_trans = run_obj["cg2:mot:det:trans"].getStatistics().mean

    sample_ap = run_obj["cg2:cs:sampleap"].getStatistics().mean

    source_ap = run_obj["cg2:cs:sourceap"].getStatistics().mean

    coll_dist = run_obj["cg2:cs:sourceaptosampleap"].getStatistics().mean

    if mag_type == "magH":
        sample_field = round(run_obj["maghfield"].getStatistics().mean, 3)
        rot = round(run_obj["p_insertion_omega"].getStatistics().mean, 3)
        dom = round(run_obj["p_insertion_omega"].getStatistics().mean, 3)
    else:
        sample_field = round(run_obj["maggfield"].getStatistics().mean, 3)
        rot = round(run_obj["srot"].getStatistics().mean, 3)
        dom = round(run_obj["CG2:Mot:MagG:rot.RBV"].getStatistics().mean, 3)

    tilt = round(run_obj["tilt"].getStatistics().mean, 3)
    sample_temp = round(run_obj["sample_temp"].getStatistics().mean, 3)

    start = run_obj["start_time"].value
    start_idx = start.find("T")
    start_time = start[:start_idx]

    start_date = start[start_idx + 1 :]

    end = run_obj["end_time"].value
    end_idx = end.find("T")
    end_time = end[:end_idx]
    end_date = end[end_idx + 1 :]

    duration = run_obj["duration"].value

    title = run_obj["run_title"].value

    monitor = run_obj["monitor"].value

    filename = os.path.basename(run_filename)

    # Regular expression to match one or more digits
    match = re.findall(r"\d+", filename)

    # If a match is found, print the numeric part
    if match:
        run_number = max(match, key=len)
        print("Numeric part:", run_number)
    else:
        print("No numeric part found in the string.")
    # print(filename)
    # filename=output_dir+'/2D/GPSANS_'+samples[idx].zfill(7)+'.dat'
    output_filename = output_dir + "/" + f"GPSANS_{run_number.zfill(7)}.dat"
    print(output_filename)
    # formatting the output file with the headers of metadata and then data.
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    f = open(output_filename, "w+")
    f.write("# " + title + "\n")
    f.write("#start_time       start_date        end_date       end_time" + "      duration    monitor\n")
    f.write(
        "# "
        + str(start_time)
        + "\t"
        + str(start_date.split(".")[0])
        + "\t"
        + str(end_date.split(".")[0])
        + "\t"
        + str(end_time)
        + "\t"
        + str(duration)
        + "\t"
        + str(monitor)
        + "\n"
    )
    f.write(
        "#lambda       dlam/lam        detector        detect_trans "
        + "       coll_dist       sample_ap       source_ap \n"
    )
    f.write(
        "# "
        + str(wavelength)
        + "\t"
        + str(wavelength_spread)
        + "\t"
        + str(detector_dist)
        + "\t"
        + str(detector_trans)
        + "\t"
        + str(coll_dist)
        + "\t"
        + str(sample_ap)
        + "\t"
        + str(source_ap)
        + "\n"
    )
    f.write("#rot   tilt   dom  sample_field   sample_temp\n")
    f.write(
        "# "
        + str(rot)
        + "\t"
        + str(tilt)
        + "\t"
        + str(dom)
        + "\t"
        + str(sample_field)
        + "\t"
        + str(sample_temp)
        + "\n"
    )  #
    f.write("#ASCII data array of 256 by 192 pixels\n\n")

    workspace = mtd[str(w_f)]
    panel_name = "detector1"
    detector_names = (
        [
            panel_name,
        ]
        if panel_name is not None
        else panel_names(workspace)
    )
    # fig = plt.figure(**figure_kwargs)
    for i_detector, detector_name in enumerate(detector_names):
        collection = TubeCollection(workspace, detector_name)
        collection = collection.sorted(view="fbfb")
        data = np.sum(np.array([tube.readY for tube in collection]), axis=-1)  # sum intensities for each pixel

    # print(data)
    # exit()
    f.write("\n".join(" ".join(map(str, x)) for x in (data)))

    #    for i in range(len(w_f.readY(0))):
    #        for j in range(w_f.getNumberHistograms()):
    #            f.write('{:.6E}\t'.format(w_f.readY(j)[i]))

    f.close()

    # mtd.clear()
    # SaveNexus(w_f, filename)

    script_end_time = time.time()
    print(script_end_time - script_start_time)
