# Modified from Volker's RC488_IPTS24666_VolkerTemplate.py

"""
Sample - /HFIR/CG3/IPTS-17240/exp318/Datafiles/BioSANS_exp318_scan0217_0001.xml (Scattering/Transmission)
Empty Beam - /HFIR/CG3/IPTS-17240/exp318/Datafiles/BioSANS_exp318_scan0220_0001.xml (For Transmission)
Beam Center - /HFIR/CG3/IPTS-17240/exp318/Datafiles/BioSANS_exp318_scan0220_0001.xml (Transmission Measurement)
Dark - /HFIR/CG3/IPTS-17240/exp318/Datafiles/BioSANS_exp318_scan0044_0001.xml (for both main and wing detectors)
"""

# IPTS and experiment
IPTS_Number = 17240
EXPERIMENT_NUMBER = 318

# Non-TimeSlice Single Configuration-
sample_identifier = ''                    # DO NOT CHANGE IF Non-TimeSlice Experiments
sample_names = ['Spice_318_217']   # DO NOT LEAVE BLANK
sample_thick = ['0.1']               # Do not repeat if sample for ALL samples
samples = [(217, 1)]       # Enter the list of runs for 'samples'
samples_trans = samples               # Enter its own list if different from 'samples' list
backgrounds = [None]                 # Do not repeat multiple times if SAME for ALL samples
backgrounds_trans = backgrounds           # Enter its own list if different from 'backgrounds' list

# Change if reducing a subset of 'samples' list
start_index = 1                     # Default start index is 1; DO NOT START FROM 'ZERO'
end_index = len(samples)          # Default is 'len(samples)'

# Setup once at the beginning of the experiment
User3LetInitial = 'whoever'                 # 3-Letter initials that identifies you in the output directory
overWrite = True                  # Option to overwrite existing data or create another folder (Default is 'False')

# ## Instrument Scientist or Local contact input below (And Expert Users)

# Advanced Settings for Data Reduction--
# Buffer clearing frequency
clearBuffer = False
refreshCycle = 25                    # loops... depending on the activities.

# Common setting for all options
scalefac = "4.05e-9"
beam_center = (220, 1)
empty_trans = (220, 1)
dark_mfname = (44, 1)
dark_wfname = (44, 1)
sens_mfname = '/HFIR/CG3/shared/Cycle488/Sens_f6368m4p0_bsSVP.nxs'
sens_wfname = '/HFIR/CG3/shared/Cycle488/Sens_f6380w1p4_bsSVP.nxs'

# Plotting range--
q_range_main = [0.003, 0.045]        # Q-range for isotropic data
q_range_wing = [0.03, 0.9]
OL_range = [0.0325, 0.0425]

# Miscellaneous settings--
base_output_directory = f'/HFIR/CG3/IPTS-{IPTS_Number}/shared/{User3LetInitial}/'
scaling_beam_radius = None
flexible_pixelsizes = True  # 'True'- if use barscan/flood information for flexible pixel sizes, else 'False'
# Make sure the barscan used sensitivity file used above if 'True'

# Plotting Options--
Plot_type = 'scalar'  # 'scalar' for isotropic and 'wedge' for anisotropic (manual or auto)
Plot_binning = 'log'  # 'log' or 'linear' Q-binning
# LINEAR BINNING
Lin1DQbins_Main = ""  # No. of bins for linear binning of 1D Main Detector, Default is 100;
# If per decade is used default is ''
Lin1DQbins_Wing = ""  # No. of bins for linear binning of 1D Wing Detector, Default is 100;
# If per decade is used default is ''
Lin2DQxy_Main = 100  # No. of bins for linear binning of 2D Main Detector, Default is 100
Lin2DQxy_Wing = 100  # No. of bins for linear binning of 2D Main Detector, Default is 100

# LOGARITHMIC BINNING
LogQbinsPerDecade_Main = 25                    # No. of bins per decade of 1D Main Detector, Default is 33
LogQbinsPerDecade_Wing = 25                    # No. of bins per decade of 1D Main Detector, Default is 33

# If time Slicing--
timeSliceExpt = False                 # 'True' if time slice experiment, else 'False'
timeSliceDuration = 60                    # Units - seconds and irrelevant if above is 'False'

# ANISOTROPIC DATA REDUCTION--
# Wedge_0...
q_range_main_wedge0 = [0.003, 0.0425]       # Q-range for anisotropic data -- wedge0
q_range_wing_wedge0 = [0.02, 0.45]
OL_range_wedge0 = [0.025, 0.04]

# Wedge_1...
q_range_main_wedge1 = [0.003, 0.0425]       # Q-range for anisotropic data -- wedge1
q_range_wing_wedge1 = [0.03, 0.45]
OL_range_wedge1 = [0.03, 0.04]

# If Manual Wedges--
wedge_min_angles = None  # If AUTOWEDGE reduction type 'None'; Irrelevant for ISOTROPIC [Wedge0_min,Wedge1_min]
wedge_max_angles = None  # If AUTOWEDGE reduction type 'None'; Irrelevant for ISOTROPIC [Wedge0_max,Wedge1_max]

# If automatic determination of wedge angles-- (To determine Wedges-TDW)
Qmin_TDW = 0.003  # Minimum Q of the Main Detector
Qmax_TDW = 0.04   # Maximum Q of the Main Detector
Qdelta_TDW = 0.01    # Q step-size (or annual ring widths); Default is 0.01. Too fine will fail autodetection.
PeakWidth_TDW = 0.5  # Wedge0 opening angle based of peak width (Default is 0.5--- 50%)
AziDelta_TDW = 1.0   # Azimuthal Angle, Phi step-size
BkgWidth_TDW = 1.0   # Wedge1 opening angle based of peak width (Default is 1.5--- 150%)
MinSigtoNoise_TDW = 2.0  # The intensity ratio between peak and background to detect a peak (Default is 2.0)

# Body of Reduction - DO NOT CHANGE....
from mantid.simpleapi import mtd  # noqa: E401
if clearBuffer:
    mtd.clear()

import numpy as np  # noqa: E401
import warnings  # noqa: E401
warnings.filterwarnings('ignore')

# get_ipython().run_line_magic('matplotlib', 'inline')
import os  # noqa: E401
import time  # noqa: E401
from drtsans.mono.spice_data import SpiceRun, map_to_nexus  # noqa: E401

# Convert SPICE scan-pt tuple to NeXus files
CG3 = 'CG3'
samples = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, samples, nexus_dir=None)
samples_trans = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, samples_trans, nexus_dir=None)
backgrounds = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, backgrounds, nexus_dir=None)
backgrounds_trans = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, backgrounds_trans, nexus_dir=None)
beam_center = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, [beam_center], nexus_dir=None)[0]
empty_trans = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, [empty_trans], nexus_dir=None)[0]
dark_mfname = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, [dark_mfname], nexus_dir=None)[0]
dark_wfname = map_to_nexus(CG3, IPTS_Number, EXPERIMENT_NUMBER, [dark_wfname], nexus_dir=None)[0]

from drtsans.mono.biosans import (load_all_files, reduce_single_configuration, plot_reduction_output,
                                  reduction_parameters, update_reduction_parameters)  # noqa: E401

# reduction parameters common to all the reduction runs to be carried out in this notebook
common_configuration = {
    "iptsNumber": IPTS_Number,
    "beamCenter": {"runNumber": beam_center},
    "emptyTransmission": {"runNumber": empty_trans},
    "configuration": {
        "outputDir": base_output_directory,
        "darkMainFileName": dark_mfname,
        "darkWingFileName": dark_wfname,
        "sensitivityMainFileName": sens_mfname,
        "sensitivityWingFileName": sens_wfname,
        "defaultMask": [{'Pixel': '1-18,239-256'}, {'Bank': '18-24,42-48'}, {'Bank': '49', 'Tube': '1'}],
        'StandardAbsoluteScale': scalefac,
        "DBScalingBeamRadius": scaling_beam_radius,
        "mmRadiusForTransmission": "",
        "absoluteScaleMethod": "standard",
        "numMainQBins": Lin1DQbins_Main,
        "numWingQBins": Lin1DQbins_Wing,
        "numMainQxQyBins": Lin2DQxy_Main,
        "numWingQxQyBins": Lin2DQxy_Wing,
        "1DQbinType": Plot_type,
        "QbinType": Plot_binning,
        "LogQBinsPerDecadeMain": LogQbinsPerDecade_Main,
        "LogQBinsPerDecadeWing": LogQbinsPerDecade_Wing,
        "useLogQBinsDecadeCenter": False,
        "useLogQBinsEvenDecade": False,
        "sampleApertureSize": 14,
        "QminMain": q_range_main[0],
        "QmaxMain": q_range_main[1],
        "QminWing": q_range_wing[0],
        "QmaxWing": q_range_wing[1],
        "overlapStitchQmin": OL_range[0],
        "overlapStitchQmax": OL_range[1],
        "usePixelCalibration": flexible_pixelsizes,
        "useTimeSlice": timeSliceExpt,
        "timeSliceInterval": timeSliceDuration,
        "WedgeMinAngles": wedge_min_angles,
        "WedgeMaxAngles": wedge_max_angles,
        "autoWedgeQmin": Qmin_TDW,
        "autoWedgeQmax": Qmax_TDW,
        "autoWedgeQdelta": Qdelta_TDW,
        "autoWedgePeakWidth": PeakWidth_TDW,
        "autoWedgeAzimuthalDelta": AziDelta_TDW,
        "autoWedgeBackgroundWidth": BkgWidth_TDW,
        "autoWedgeSignalToNoiseMin": MinSigtoNoise_TDW,
        "wedge1QminMain": q_range_main_wedge0[0],
        "wedge1QmaxMain": q_range_main_wedge0[1],
        "wedge1QminWing": q_range_wing_wedge0[0],
        "wedge1QmaxWing": q_range_wing_wedge0[1],
        "wedge1overlapStitchQmin": OL_range_wedge0[0],
        "wedge1overlapStitchQmax": OL_range_wedge0[1],
        "wedge2QminMain": q_range_main_wedge1[0],
        "wedge2QmaxMain": q_range_main_wedge1[1],
        "wedge2QminWing": q_range_wing_wedge1[0],
        "wedge2QmaxWing": q_range_wing_wedge1[1],
        "wedge2overlapStitchQmin": OL_range_wedge1[0],
        "wedge2overlapStitchQmax": OL_range_wedge1[1],
    }
}

common_configuration_full = reduction_parameters(common_configuration, 'BIOSANS', validate=False)
# pretty_print(common_configuration_full)

if len(backgrounds) == 1 and len(samples) > len(backgrounds):
    backgrounds = backgrounds*len(samples)
if len(backgrounds_trans) == 1 and len(samples_trans) > len(backgrounds_trans):
    backgrounds_trans = backgrounds_trans*len(samples_trans)
if len(sample_thick) == 1 and len(samples) > len(sample_thick):
    sample_thick = sample_thick*len(samples)

# Checking if output directory exists, if it doesn't, creates the folder
# Also, if do not overwrite, then makes sure the directory does not exists.
output_dir = base_output_directory
if not overWrite:
    suffix = 0
    while os.path.exists(output_dir):
        output_dir = base_output_directory[0, len(base_output_directory)-2] + "_" + str(suffix) + "/"
        suffix += 1

if not timeSliceExpt and sample_identifier is not '':
    if sample_identifier is not '':
        output_dir = base_output_directory+str(sample_identifier)+"/"
        change_outputdir = {
            'configuration': {
                'outputDir': output_dir,
            },
        }
        common_configuration_full = update_reduction_parameters(common_configuration_full,
                                                                change_outputdir,
                                                                validate=False)
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)


start_time = time.time()
# Loop for samples
for i in range(start_index-1, end_index):
    start_time_loop = time.time()
    if timeSliceExpt:
        output_dir = base_output_directory+"timeslice/t"+str(timeSliceDuration)+"/"+sample_names[i]+"/"
        timeslice_outputdir = {
            'configuration': {
                'outputDir': output_dir,
            },
        }
        common_configuration_full = update_reduction_parameters(common_configuration_full,
                                                                timeslice_outputdir,
                                                                validate=False)
        for subfolder in ['1D', '2D']:
            output_folder = os.path.join(output_dir, subfolder)
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)

    print("Reducing...", samples[i], ":", sample_names[i], '\n')
    run_data = {
        'sample': {
            'runNumber': samples[i],
            'thickness': sample_thick[i],
            'transmission': {'runNumber': samples_trans[i]}
        },
        'background': {
            'runNumber': backgrounds[i],
            'transmission': {'runNumber': backgrounds_trans[i]}
        },
        'outputFileName': f'r{samples[i]}_{sample_names[i]}',
    }

    # Update our common settings with the particulars of the current reduction
    reduction_input = update_reduction_parameters(common_configuration_full, run_data, validate=True)
    # pretty_print(reduction_input)
    reduction_input['configuration']['WedgeMinAngles'] = wedge_min_angles
    reduction_input['configuration']['WedgeMaxAngles'] = wedge_max_angles
    loaded = load_all_files(reduction_input)
    out = reduce_single_configuration(loaded, reduction_input)
    plot_reduction_output(out, reduction_input)

    print('\nloop_'+str(i+1)+": ", time.time()-start_time_loop)

    if np.remainder(i, refreshCycle) == 0 and i > 0:
        mtd.clear()

print('Total Time : ', time.time()-start_time)

from mantid import mtd  # noqa: E401
mtd.getObjectNames()
