from __future__ import print_function

import numpy as np
from scipy import stats

from ornl.sans.momentum_transfer \
    import MomentumTransfer as MomentumTransferMain


# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')


class MomentumTransfer(MomentumTransferMain):

    def __iadd__(self, other):
        """This is an overload for `+=` operator.
        Usefull for EQSANS when calculating I(qx,qy) or I(q) from several bins.

        Parameters
        ----------
        other : MomentumTransfer class

        """
        self.qx = np.concatenate(self.qx, other.qx)
        self.qy = np.concatenate(self.qy, other.qy)
        self.dqx = np.concatenate(self.dqx, other.dqx)
        self.dqy = np.concatenate(self.dqy, other.dqy)
        self.i = np.concatenate(self.i, other.i)
        self.i_sigma = np.concatenate(self.i_sigma, other.i_sigma)