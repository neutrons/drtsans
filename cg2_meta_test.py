# From:
#   gpsans_reduction_1config.ipynb

#USER Input here with scan numbers etc.
samples = ['9166', '9167', '9176']
samples_trans = ['9178', '9179', '9188']
sample_thick = ['0.1']*3
bkgd = ['9165','9165','9165']
bkgd_trans = ['9177', '9177', '9177']


sample_names = ["Al4", "PorasilC3", "PTMA-15"]
##########################Local contact adds json file location####################
json_file = '/HFIR/CG2/shared/UserAcceptance/LDS_metadata/gpsans_reduction.json'

###############Reduction start here, no user input needed below ########################

import os
import json

# chekcing if output directory exists, if it doesn't, creates the folder


with open(json_file) as f:
    reduction_input = json.load(f)

output_dir = '/SNS/users/wzz/Projects/SANS/sans-backend/temp/'
reduction_input["configuration"]["outputDir"] = output_dir
for subfolder in ['1D', '2D']:
    output_folder = os.path.join(output_dir, subfolder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

from drtsans.mono.gpsans import load_all_files, reduce_single_configuration, plot_reduction_output
import time

# About SWD and SDD
swd = reduction_input["configuration"]["SampleToSi"]
sdd = reduction_input["configuration"]["SampleDetectorDistance"]
print('JSON: SWD = {}, SDD = {}'.format(swd, sdd))

start_time = time.time()
for i in range(len(samples)):
    reduction_input["runNumber"] = samples[i]
    reduction_input["transmission"]["runNumber"] = samples_trans[i]
    reduction_input["background"]["runNumber"] = bkgd[i]
    reduction_input["background"]["transmission"]["runNumber"] = bkgd_trans[i]
    reduction_input["outputFilename"] = sample_names[i]
    reduction_input["thickness"] = sample_thick[i]
    loaded = load_all_files(reduction_input)
    out = reduce_single_configuration(loaded, reduction_input)
    plot_reduction_output(out, reduction_input, loglog=False)

end_time = time.time()
print(end_time - start_time)



