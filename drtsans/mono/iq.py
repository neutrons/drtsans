import numpy as np

from drtsans.mono import momentum_transfer
from drtsans.iq import \
    MomentumTransfer as MomentumTransferMain


# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')


class MomentumTransfer(MomentumTransferMain):
    def __init__(self,
                 input_workspace=None,
                 component_name="detector1",
                 out_ws_prefix="ws"):
        super(MomentumTransfer, self).__init__(momentum_transfer, input_workspace,
                                               component_name=component_name,
                                               out_ws_prefix=out_ws_prefix)
