
from ornl.sans.sns.eqsans.parameters import get_parameters

from mantid.simpleapi import Load


# beam center file
FILENAME = os.path.join(os.getenv('DATA_DIRECTORY'), 'eqsans', 'EQSANS_68183_event.nxs'),

def direct_beam_center():
    params = get_parameters(68183)
    data = Load(FILENAME)
    

