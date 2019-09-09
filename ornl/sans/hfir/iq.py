from __future__ import print_function

import numpy as np

from ornl.sans.hfir import resolution
from ornl.sans.iq import \
    MomentumTransfer as MomentumTransferMain


# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')


class MomentumTransfer(MomentumTransferMain):
    def __init__(self,
                 input_workspace=None,
                 component_name="detector1",
                 out_ws_prefix="ws"):
        super(MomentumTransfer, self).__init__(resolution, input_workspace,
                                               component_name=component_name,
                                               out_ws_prefix=out_ws_prefix)
