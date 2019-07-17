from __future__ import print_function

import numpy as np

from ornl.sans.momentum_transfer \
    import MomentumTransfer as MomentumTransferMain


class MomentumTransfer(MomentumTransferMain):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __iadd__(self, other):
        """This is an overload for `+=` operator.
        Usefull for EQSANS when calculating I(qx,qy) or I(q) from several bins.

        Parameters
        ----------
        other : MomentumTransfer class

        """
        self.qx = np.concatenate((self.qx, other.qx))
        self.qy = np.concatenate((self.qy, other.qy))
        self.dqx = np.concatenate((self.dqx, other.dqx))
        self.dqy = np.concatenate((self.dqy, other.dqy))
        self.i = np.concatenate((self.i, other.i))
        self.i_sigma = np.concatenate((self.i_sigma, other.i_sigma))

        return self
