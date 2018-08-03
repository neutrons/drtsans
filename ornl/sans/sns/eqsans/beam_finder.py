
import os
import sys

from ornl.sans.sns.eqsans.parameters import get_parameters


from pprint import pprint

from dotenv import load_dotenv
load_dotenv()

sys.path.append(os.getenv("MANTID_PATH"))
from mantid.simpleapi import (
    Load, SANSMaskDTP, Integration, FindCenterOfMassPosition)

# beam center file
FILENAME = os.path.abspath(os.path.join(os.getenv('DATA_DIRECTORY'), 'eqsans', 'EQSANS_68183_event.nxs'))
TUBES_TO_MASK = "1,48,53,54,85,123,130,137" # From 1 to 192
    

def direct_beam_center():
    ws = Load(FILENAME)
    SANSMaskDTP(InputWorkspace=ws, Tube=TUBES_TO_MASK)
    # Flatten TOF
    ws_flattened = Integration(InputWorkspace=ws)
    center = FindCenterOfMassPosition(InputWorkspace=ws_flattened)
    center_x, center_y = center
    return center_x, center_y
    # TO Center the instrument
    #MoveInstrumentComponent(Workspace=ws, ComponentName='detector1', X=-center_x, Y=-center_y)




if __name__ == "__main__":
    direct_beam_center()