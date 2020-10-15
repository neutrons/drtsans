# Convert BIOSANS SPICE file to event NeXus
# All the output file is supposed to be written to
# /HFIR/CG3/IPTS-{ipts_number}/shared/spice_nexus/CG3_{exp}{scan}{pt}'
# and make up a unique run number from experiment number, scan number and pt number


# Set SPICE files information
# The following example is for sensitivities preparation

ipts = 17241
exp = 549
scan_pt_list = zip([9, 10, 20, 16, 22], [1] * 5)

# ----------------------------------------------------------------------------------
# TRY NOT TO TOUCH THIS PART
# ----------------------------------------------------------------------------------
# Template event nexus file for instrument geometry and etc
TEMPLATE_EVENT_NEXUS = "/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/biosans/CG3_5705.nxs.h5"

# ----------------------------------------------------------------------------------
# DON'T TOUCH ANYTHING BELOW THIS LINE
# ----------------------------------------------------------------------------------
from drtsans.mono.biosans.cg3_spice_to_nexus import convert_spice_to_nexus  # noqa: E401

nexus_names = set()
for scan_num, pt_num in scan_pt_list:
    nexus = convert_spice_to_nexus(ipts, exp, scan_num, pt_num, TEMPLATE_EVENT_NEXUS,
                                   masked_detector_pixels=[70911],
                                   output_dir=f'/HFIR/CG3/IPTS-{ipts}/shared/Exp{exp}')
    nexus_names.add(nexus)

# Check
from mantid.simpleapi import LoadEventNexus  # noqa: E401
print(f'\n\n-----  Verification -------\n\n')
for nexus_name in nexus_names:
    try:
        print(f'Loading {nexus_name}')
        ws = LoadEventNexus(Filename=nexus_name, LoadNexusInstrumentXML=True)
        print(f'{nexus_name}: {ws.getNumberHistograms()}, {ws.getNumberEvents()}')
    except RuntimeError as run_err:
        print(f'Failed to laod {nexus_name} due to {run_err}')
