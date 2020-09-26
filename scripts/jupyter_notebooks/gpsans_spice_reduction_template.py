# This notebook is a template for reducing GPSANS legacy data with drtsans
# It is a 1-sample 1-configuration reduction

# USER INPUT
ipts_number = 828
exp_number = 280

# Note:
# Default directory for converted SPICE data files is /HFIR/CG2/IPTS-828/shared/Exp280/
# Here is the convention for a SPICE with an example
# Exp = 280, Scan = 5, Pt = 29: CG2_028000050029.nxs.h5

output_directory = f'/HFIR/CG2/IPTS-{ipts_number}/shared/reduced/Exp{exp_number}/'

samples = [(35, 1)]
samples_trans = [(27, 1)]
sample_thick = ['0.1']
sample_names = ['Porasil B']
bkgd = [(34, 1)]
bkgd_trans = [(26, 1)]

empty_trans = [(28, 1)]
beam_center = [(20, 1)]

# q range to use to clean 1D curve of each configuration
q_range = [None, None]

# STAFF INPUT
use_mask_file = True
mask_file_name = f'/HFIR/CG2/IPTS-{ipts_number}/shared/pixel_calibration/mask_pixel_map.nxs'
use_dark_file = False
dark_file_name = ""
block_beam = (9, 1)
use_mask_back_tubes = False
wavelength = None
wavelength_spread = None
wedge_min_angles = None
wedge_max_angles = None


# Convert from SPICE naming system to converted Nexus files
# DO NOT TOUCH
def map_to_nexus(ipts, exp, scan_pt_list):
    import os
    nexus_list = list()
    for scan, pt in scan_pt_list:
        nexus_file = f'/HFIR/CG2/IPTS-{ipts}/shared/Exp{exp}/CG2_{exp:04}{scan:04}{pt:04}.nxs.h5'
        if os.path.exists(nexus_file) is False:
            print(f'WARNING: Scan {scan} Pt {pt}: event Nexus {nexus_file} does not exist.')
        nexus_list.append(nexus_file)
    return nexus_list


samples = map_to_nexus(ipts_number, exp_number, samples)
samples_trans = map_to_nexus(ipts_number, exp_number, samples_trans)
bkgd = map_to_nexus(ipts_number, exp_number, bkgd)
bkgd_trans = map_to_nexus(ipts_number, exp_number, bkgd_trans)
empty_trans = map_to_nexus(ipts_number, exp_number, empty_trans)
beam_center = map_to_nexus(ipts_number, exp_number, beam_center)
block_beam = map_to_nexus(ipts_number, exp_number, [block_beam])[0]
# --- END OF NO TOUCH ZONE ---

see_full_json = True  # To display full configuration json (True/False)
use_log_2d_binning = False
use_log_1d = True
common_configuration = {
    "iptsNumber": ipts_number,
    "emptyTransmission": {"runNumber": empty_trans},
    "beamCenter": {"runNumber": beam_center},
    "configuration": {
        "outputDir": output_directory,
        "darkFileName": dark_file_name,
        "sensitivityFileName": '/HFIR/CG2/shared/drt_sensitivity/sens_fc488_bar.nxs',
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
        "blockedBeamRunNumber": block_beam,
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

# Never touch!  drtsans specific
import warnings
warnings.filterwarnings('ignore')
import os
import json
# jupyter only:from pprint import pprint as pretty_print
import time
from drtsans.mono.gpsans import (load_all_files, reduce_single_configuration, plot_reduction_output,
                                 reduction_parameters, update_reduction_parameters)
from matplotlib.colors import LogNorm
# jupyter only:% matplotlib
# jupyter only:inline

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
            "blockedBeamRunNumber": block_beam,
            "maskFileName": mask_file_name,
            "darkFileName": dark_file_name,
        }
    }

    # Update our common settings with the particulars of the current reduction
    reduction_input = update_reduction_parameters(common_configuration_full, run_data, validate=True)

    # Begin redution. Be sure to validate the parameters before.
    # We validated the parameters when we updated them. Otherwise you can invoke the validation like this:
    # reduction_input = validate_reduction_parameters(reduction_input)

    loaded = load_all_files(reduction_input, path=f'/HFIR/CG2/IPTS-{ipts_number}/shared/Exp{exp_number}')
    out = reduce_single_configuration(loaded, reduction_input)
    plot_reduction_output(out, reduction_input, loglog=use_log_1d, imshow_kwargs=log_flag)

    # Save the reduction parameters of each reduction session to a JSON file
    output_dir = reduction_input['configuration']['outputDir']
    output_json_file = os.path.join(output_dir, f'{sample_names[i]}.json')  # full path to the JSON file
    with open(output_json_file, 'w') as file_handle:
        json.dump(reduction_input, file_handle, indent=2)

end_time = time.time()
print("Total run time: {}s".format(end_time - start_time))