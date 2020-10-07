# This notebook is a template for reducing GPSANS legacy data with drtsans
# It is a 1-sample 1-configuration reduction
# ----------------------  HEADER: NO TOUCH  ------------------------------
from drtsans.mono.spice_data import map_to_nexus
# Never touch!  drtsans specific
import warnings  # noqa E402
warnings.filterwarnings('ignore')
import os  # noqa E402
import json  # noqa E402
# jupyter only:from pprint import pprint as pretty_print
import time # noqa E402
from drtsans.mono.gpsans import (load_all_files, reduce_single_configuration, plot_reduction_output,
                                 reduction_parameters, update_reduction_parameters) # noqa E402
from matplotlib.colors import LogNorm # noqa E402
# jupyter only:% matplotlib
# jupyter only:inline
# ----------------------  HEADER: END      ------------------------------

# USER INPUT
CG2 = 'CG2'
ipts_number = 828
exp_number = 280

# Default directory for converted SPICE data files is /HFIR/CG2/IPTS-828/shared/Exp280/
# Here is the convention for a SPICE with an example
# Exp = 280, Scan = 5, Pt = 29: CG2_028000050029.nxs.h5
output_directory = f'/HFIR/{CG2}/IPTS-{ipts_number}/shared/reduced/Exp{exp_number}/'

# Sample
samples = [(35, 1)]
samples_trans = [(27, 1)]
sample_thick = ['0.1']
sample_names = ['Porasil_B']
# Background
bkgd = [(34, 1)]
bkgd_trans = [(26, 1)]
# Transmission/empty beam
empty_trans = [(28, 1)]
# Beam center
beam_center = [(20, 1)]

# q range to use to clean 1D curve of each configuration
q_range = [None, None]

# STAFF INPUT
# Block beam
block_beam = (9, 1)

# Others: mask, dark file, ...
use_mask_file = True
mask_file_name = f'/HFIR/{CG2}/IPTS-{ipts_number}/shared/pixel_calibration/mask_pixel_map.nxs'
use_dark_file = False
dark_file_name = ""
use_mask_back_tubes = False

# Wave length
wavelength = None
wavelength_spread = None

# Wedge
wedge_min_angles = None
wedge_max_angles = None

sensitivity_file = '/HFIR/CG2/shared/drt_sensitivity/sens_CG2_spice_bar.nxs'
see_full_json = True  # To display full configuration json (True/False)
use_log_2d_binning = False
use_log_1d = True
common_configuration = {
    "iptsNumber": ipts_number,
    "emptyTransmission": {"runNumber": map_to_nexus(CG2, ipts_number, exp_number, empty_trans)},
    "beamCenter": {"runNumber": map_to_nexus(CG2, ipts_number, exp_number, beam_center)},
    "configuration": {
        "outputDir": output_directory,
        "darkFileName": dark_file_name,
        "sensitivityFileName": sensitivity_file,
        "DBScalingBeamRadius": 40,
        "sampleApertureSize": 8,
        "mmRadiusForTransmission": 40,
        "absoluteScaleMethod": "direct_beam",
        "numQxQyBins": 256,
        "1DQbinType": "scalar",
        "QbinType": "log",
        "numQBins": "",
        "LogQBinsPerDecade": 33,
        "useLogQBinsDecadeCenter": True,
        "useLogQBinsEvenDecade": False,
        "wavelength": wavelength,
        "wavelengthSpread": wavelength_spread,
        "maskFileName": mask_file_name,
        'WedgeMinAngles': wedge_min_angles,
        'WedgeMaxAngles': wedge_max_angles,
        "AnnularAngleBin": 2.0,
        "Qmin": 0.0028,
        "Qmax": 0.0035,
        "useSubpixels": True,
        "subpixelsX": 5,
        "subpixelsY": 5,
        "useTimeSlice": False,
        "useLogSlice": False,
        "logSliceName": "",
        "logSliceInterval": '',
    }
}

# ------------------- Never touch!  drtsans specific -----------------------------

# convert SPICE to Nexus
samples = map_to_nexus(CG2, ipts_number, exp_number, samples)
samples_trans = map_to_nexus(CG2, ipts_number, exp_number, samples_trans)
bkgd = map_to_nexus(CG2, ipts_number, exp_number, bkgd)
bkgd_trans = map_to_nexus(CG2, ipts_number, exp_number, bkgd_trans)

if not use_mask_file:
    mask_file_name = ""
if not use_dark_file:
    dark_file_name = ""

if use_log_2d_binning:
    log_flag = {"norm": LogNorm()}
else:
    log_flag = {'vmin': 0, 'vmax': 100}

# Add on the other reduction parameters with their default values (most will be empty)
common_configuration_full = reduction_parameters(common_configuration, 'GPSANS', validate=False)

# Create output directory

output_dir = common_configuration_full['configuration']['outputDir']
for subfolder in ['1D', '2D']:
    output_folder = os.path.join(output_dir, subfolder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

start_time = time.time()
for i in range(len(samples)):
    # Settings particular to each reduction session
    run_data = {
        'sample': {
            'runNumber': samples[i],
            'thickness': sample_thick[i],
            'transmission': {'runNumber': samples_trans[i]}
        },
        'background': {
            'runNumber': bkgd[i],
            'transmission': {'runNumber': bkgd_trans[i]}
        },
        'outputFileName': sample_names[i],
        'configuration': {
            "Qmin": q_range[0],
            "Qmax": q_range[1],
            "useMaskBackTubes": use_mask_back_tubes,
            "blockedBeamRunNumber": map_to_nexus(CG2, ipts_number, exp_number, [block_beam])[0],
            "maskFileName": mask_file_name,
            "darkFileName": dark_file_name,
        }
    }

    # Update our common settings with the particulars of the current reduction
    reduction_input = update_reduction_parameters(common_configuration_full, run_data, validate=True)

    # Begin redution. Be sure to validate the parameters before.
    # We validated the parameters when we updated them. Otherwise you can invoke the validation like this:
    # reduction_input = validate_reduction_parameters(reduction_input)

    loaded = load_all_files(reduction_input, path=f'/HFIR/{CG2}/IPTS-{ipts_number}/shared/Exp{exp_number}')
    out = reduce_single_configuration(loaded, reduction_input)
    plot_reduction_output(out, reduction_input, loglog=use_log_1d, imshow_kwargs=log_flag)

    # Save the reduction parameters of each reduction session to a JSON file
    output_dir = reduction_input['configuration']['outputDir']
    output_json_file = os.path.join(output_dir, f'{sample_names[i]}.json')  # full path to the JSON file
    with open(output_json_file, 'w') as file_handle:
        json.dump(reduction_input, file_handle, indent=2)

end_time = time.time()
print("Total run time: {}s".format(end_time - start_time))
