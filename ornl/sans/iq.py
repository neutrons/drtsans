from __future__ import print_function

from ast import literal_eval
from string import Template

import numpy as np
from scipy import stats

import mantid
from mantid.kernel import logger
from mantid.simpleapi import CreateEmptyTableWorkspace, CreateWorkspace
from ornl.sans.detector import Component


# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')


class MomentumTransfer:
    '''
    Momentum Transfer class
    Olds arrays for qx, qy, dqx, dqy, i, i_sigma
    q1d and q2d are calculated from there arrays
    '''

    # For now the detector dims in Table Workspace will be in the comment
    # as below (there is no other form of adding key=value to table WS)
    DETECTOR_DIMENSIONS_TEMPLATE = "detector_dimensions=($dim_x,$dim_y)"

    qx, qy, dqx, dqy, i, i_sigma = np.empty(0), np.empty(0), np.empty(0), \
        np.empty(0), np.empty(0), np.empty(0)
    prefix = None
    component = None
    detector_dims = None

    def __init__(self, resolution, input_workspace=None,
                 component_name="detector1", out_ws_prefix="ws"):

        self.resolution = resolution
        self.prefix = out_ws_prefix

        if isinstance(input_workspace, mantid.dataobjects.TableWorkspace):
            self._load_table_workspace(input_workspace)

        elif isinstance(input_workspace, mantid.dataobjects.Workspace2D):
            self.component = Component(input_workspace, component_name)
            self.detector_dims = (self.component.dim_x, self.component.dim_y)
            self._initialize_qs(input_workspace)

    def _load_table_workspace(self, input_workspace):
        data = input_workspace.toDict()
        self.qx = np.array(data['Qx'])
        self.qy = np.array(data['Qy'])
        self.dqx = np.array(data['dQx'])
        self.dqy = np.array(data['dQy'])
        self.i = np.array(data["I"])
        self.i_sigma = np.array(data["Sigma(I)"])

        # Gets the detector dimensions: "detector_dimensions=($dim_x,$dim_y)"
        if self.DETECTOR_DIMENSIONS_TEMPLATE.split("=")[0] != \
                input_workspace.getComment().split("=")[0]:
            logger.error("Can not get the detector dimensions!")
        else:
            self.detector_dims = literal_eval(
                input_workspace.getComment().split("=")[1])

    def _remove_monitors(self, data):
        return data[self.component.first_index:self.component.first_index +
                    self.component.dims]

    def _mask_pixels(self, data, masked_pixels):
        return np.ma.MaskedArray(data,
                                 masked_pixels,
                                 dtype=np.float,
                                 fill_value=np.nan)

    def _initialize_qs(self, input_workspace):

        assert input_workspace.blocksize() == 1  # sanity check: only 1 bin

        masked_pixels = self.component.masked_ws_indices()

        # 1D arrays
        self.qx, self.qy, self.dqx, self.dqy = \
            self.resolution.q_resolution_per_pixel(input_workspace)
        self.i = input_workspace.extractY().ravel()
        self.i_sigma = input_workspace.extractE().ravel()
        # Ravel just in case! For EQSANS for example is needed!
        self.qx, self.qy, self.dqx, self.dqy = \
            self.qx.ravel(), self.qy.ravel(), self.dqx.ravel(), \
            self.dqy.ravel()
        # Get rid of the monitors; from 49154 to 49152 spectra
        self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma = [
            self._remove_monitors(d) for d in
            [self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma]
        ]

        # Create numpy mask arrays with the masked pixels
        self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma = [
            self._mask_pixels(d, masked_pixels) for d in
            [self.qx, self.qy, self.dqx, self.dqy, self.i, self.i_sigma]
        ]

    def _create_table_ws(self, suffix):
        table_iq = CreateEmptyTableWorkspace(OutputWorkspace=self.prefix +
                                             suffix)
        table_iq.addColumn(type="float", name="Qx")
        table_iq.addColumn(type="float", name="Qy")
        table_iq.addColumn(type="float", name="dQx")
        table_iq.addColumn(type="float", name="dQy")
        table_iq.addColumn(type="float", name="I")
        table_iq.addColumn(type="float", name="Sigma(I)")

        for qx_i, qy_i, dqx_i, dqy_i, i_i, i_sigma_i in zip(
                self.qx.tolist(), self.qy.tolist(), self.dqx.tolist(),
                self.dqy.tolist(), self.i.tolist(), self.i_sigma.tolist()):
            nextRow = {
                'Qx': qx_i,
                'Qy': qy_i,
                'dQx': dqx_i,
                'dQy': dqy_i,
                "I": i_i,
                "Sigma(I)": i_sigma_i,
            }
            table_iq.addRow(nextRow)

        template = Template(self.DETECTOR_DIMENSIONS_TEMPLATE)
        table_iq.setComment(
            template.substitute(dim_x=self.detector_dims[0],
                                dim_y=self.detector_dims[1]))

        return table_iq

    def q2d(self, suffix="_iqxqy_table"):
        return self._create_table_ws(suffix)

    def bin_into_q2d(self, bins=None, suffix="_iqxqy"):
        """Bin the data into Q 2D for visualization only!
        No uncertainties data!!!

        Parameters
        ----------

        bins : int or array_like or [int, int] or [array, array], optional
            The bin specification:
                - If int, the number of bins for the two dimensions
                (nx=ny=bins).
                - If array_like, the bin edges for the two dimensions
                (x_edges=y_edges=bins).
                - If [int, int], the number of bins in each dimension
                (nx, ny = bins).
                - If [array, array], the bin edges in each dimension
                (x_edges, y_edges = bins).
                - A combination [int, array] or [array, int], where int is the
                number of bins and array is the bin edges.
            If None the detector dimensions are given:
            bins = [det.dim_x, det.dim_y]

        Returns
        -------
        string, Workspace2D

        """

        if bins is None:
            bins = [self.component.dim_x, self.component.dim_y]

        # Number of bins in Qx Qy is the number of pixels in X and Y
        counts_qx_qy_weights, qx_bin_edges, qy_bin_edges = np.histogram2d(
            self.qx, self.qy, bins=bins)
        counts_qx_qy, qx_bin_edges, qy_bin_edges = np.histogram2d(
            self.qx, self.qy, bins=bins, weights=self.i)
        counts_qx_qy /= counts_qx_qy_weights

        qy_bin_centers = (qy_bin_edges[1:] + qy_bin_edges[:-1]) / 2.0

        # When doing histogram2d the masks are gone
        iqxqy_ws = CreateWorkspace(DataX=np.tile(qx_bin_edges,
                                                 len(qy_bin_centers)),
                                   DataY=np.fliplr(counts_qx_qy.T),
                                   NSpec=len(qy_bin_centers),
                                   UnitX='MomentumTransfer',
                                   VerticalAxisUnit='MomentumTransfer',
                                   VerticalAxisValues=qy_bin_centers,
                                   OutputWorkspace=self.prefix + suffix)

        return iqxqy_ws.name(), iqxqy_ws

    @staticmethod
    def _bin_into_q1d(q, dq, i, i_sigma, prefix, bins=100,
                      statistic='mean', suffix="_iq"):

        intensity_statistic, q_bin_edges, _ = stats.binned_statistic(
            q, i, statistic=statistic, bins=bins)

        sigma_statistic, q_bin_edges, _ = stats.binned_statistic(
            q,
            i_sigma,
            statistic=lambda array_1d: np.sqrt(np.sum(np.square(array_1d))
                                               ) / len(array_1d),
            bins=bins)

        _, dq_bin_edges, _ = \
            stats.binned_statistic(dq, i, statistic=statistic, bins=bins)

        dq_bin_centers = (dq_bin_edges[1:] + dq_bin_edges[:-1]) / 2.0

        iq = CreateWorkspace(
            DataX=np.array([q_bin_edges]),
            DataY=np.array([intensity_statistic]),
            DataE=np.array([sigma_statistic]),
            Dx=dq_bin_centers,  # bin centers!!
            NSpec=1,
            UnitX='MomentumTransfer',
            YUnitLabel='Counts',
            OutputWorkspace=prefix + suffix)

        return iq.name(), iq

    def bin_into_q1d(self, bins=100, statistic='mean', suffix="_iq"):
        """ Calculates: I(Q) and Dq
        The ws_* input parameters are the output workspaces from bin_into_q2d

        Parameters
        ----------

        bins : int or sequence of scalars, optional
            See `scipy.stats.binned_statistic`.
            If `bins` is an int, it defines the number of equal-width bins in
            the given range (10 by default).  If `bins` is a sequence, it
            defines the bin edges, including the rightmost edge, allowing for
            non-uniform bin widths.  Values in `x` that are smaller than lowest
            bin edge areassigned to bin number 0, values beyond the highest bin
            are assigned to ``bins[-1]``.  If the bin edges are specified,
            the number of bins will be, (nx = len(bins)-1).
        statistic : str, optional
            See `scipy.stats.binned_statistic`.
            The statistic to compute, by default 'mean'
            The following statistics are available:
            * 'mean' : compute the mean of values for points within each bin.
                Empty bins will be represented by NaN.
            * 'std' : compute the standard deviation within each bin. This
                is implicitly calculated with ddof=0.
            * 'median' : compute the median of values for points within each
                bin. Empty bins will be represented by NaN.
            * 'count' : compute the count of points within each bin.  This is
                identical to an unweighted histogram.  `values` array is not
                referenced.
            * 'sum' : compute the sum of values for points within each bin.
                This is identical to a weighted histogram.
            * 'min' : compute the minimum of values for points within each bin.
                Empty bins will be represented by NaN.
            * 'max' : compute the maximum of values for point within each bin.
                Empty bins will be represented by NaN.
            * function : a user-defined function which takes a 1D array of
                values, and outputs a single numerical statistic. This function
                will be called on the values in each bin.  Empty bins will be
                represented by function([]), or NaN if this returns an error.
        out_ws_prefix : str, optional
            The prefix of the workspace created in Mantid, by default "ws"

        Returns
        -------
        tuple
            tuple (workspace name, workspace)

        """

        q = np.sqrt(np.square(self.qx) + np.square(self.qy))
        dq = np.sqrt(np.square(self.dqx) + np.square(self.dqy))

        return MomentumTransfer._bin_into_q1d(q, dq, self.i, self.i_sigma,
                                              self.prefix, bins, statistic,
                                              suffix)

    def bin_wedge_into_q1d(self, phi_0=0, phi_aperture=30, bins=100,
                           statistic='mean', suffix="_wedge_iq"):
        """
        Wedge calculation and integration

        Calculates: I(Q) and Dq
        The ws_* input parameters are the output workspaces from bin_into_q2d

        Parameters
        ----------
        phi_0 : int, optional
            Where to start the wedge, by default 0
        phi_aperture : int, optional
            Aperture of the wedge, by default 30
        bins : int or sequence of scalars, optional
            See `scipy.stats.binned_statistic`.
            If `bins` is an int, it defines the number of equal-width bins in
            the given range (10 by default).  If `bins` is a sequence, it
            defines the bin edges, including the rightmost edge, allowing for
            non-uniform bin widths.  Values in `x` that are smaller than lowest
            bin edge areassigned to bin number 0, values beyond the highest bin
            are assigned to ``bins[-1]``.  If the bin edges are specified,
            the number of bins will be, (nx = len(bins)-1).
        statistic : str, optional
            See `scipy.stats.binned_statistic`.
            The statistic to compute, by default 'mean'
            The following statistics are available:
            * 'mean' : compute the mean of values for points within each bin.
                Empty bins will be represented by NaN.
            * 'std' : compute the standard deviation within each bin. This
                is implicitly calculated with ddof=0.
            * 'median' : compute the median of values for points within each
                bin. Empty bins will be represented by NaN.
            * 'count' : compute the count of points within each bin.  This is
                identical to an unweighted histogram.  `values` array is not
                referenced.
            * 'sum' : compute the sum of values for points within each bin.
                This is identical to a weighted histogram.
            * 'min' : compute the minimum of values for points within each bin.
                Empty bins will be represented by NaN.
            * 'max' : compute the maximum of values for point within each bin.
                Empty bins will be represented by NaN.
            * function : a user-defined function which takes a 1D array of
                values, and outputs a single numerical statistic. This function
                will be called on the values in each bin.  Empty bins will be
                represented by function([]), or NaN if this returns an error.
        out_ws_prefix : str, optional
            The prefix of the workspace created in Mantid, by default "ws"

        Returns
        -------
        tuple
            tuple (workspace name, workspace)

        """

        angle_rad = np.arctan2(self.qy, self.qx)

        # plt.scatter(mt.qx, mt.qy, c=angle_rad)
        # plt.colorbar()
        # plt.show()

        # # This is just to show the angle
        # CreateWorkspace(
        #     DataX=self.qx,
        #     DataY=angle_grid,
        #     DataE=np.sqrt(angle_grid),
        #     NSpec=256,
        #     UnitX='MomentumTransfer',
        #     VerticalAxisUnit='MomentumTransfer',
        #     VerticalAxisValues=self.qy,
        #     OutputWorkspace=out_ws_prefix+"_angle",
        # )

        # Let's work in radians
        phi_0_rad = np.deg2rad(phi_0)
        phi_aperture_rad = np.deg2rad(phi_aperture)

        phi_aperture_min_rad = phi_0_rad - phi_aperture_rad / 2
        phi_aperture_max_rad = phi_0_rad + phi_aperture_rad / 2
        # opposite 180 degrees apart
        phi_aperture_min_pi_rad = phi_aperture_min_rad + np.pi
        phi_aperture_max_pi_rad = phi_aperture_max_rad + np.pi

        condition1 = (angle_rad > phi_aperture_min_rad) & \
            (angle_rad < phi_aperture_max_rad)
        # make angle > np.pi varying between np.pi and 2*np.pi, rather than the
        # initial -np.pi to np.pi
        angle_rad[angle_rad < 0] = 2 * np.pi + angle_rad[angle_rad < 0]
        condition2 = (angle_rad > phi_aperture_min_pi_rad) & \
            (angle_rad < phi_aperture_max_pi_rad)
        # True where the wedge is located, otherwise false
        condition = condition1 | condition2

        # # This is just to show the condition
        # CreateWorkspace(
        #     DataX=qx_bin_edges_grid,
        #     DataY=condition,
        #     DataE=np.sqrt(angle_grid),
        #     NSpec=256,
        #     UnitX='MomentumTransfer',
        #     VerticalAxisUnit='MomentumTransfer',
        #     VerticalAxisValues=qy_bin_centers,
        #     OutputWorkspace=out_ws_prefix+"_condition",
        # )

        q = np.sqrt(np.square(self.qx) + np.square(self.qy))
        q = q[condition]

        i = self.i[condition]
        i_sigma = self.i_sigma[condition]

        dq = np.sqrt(np.square(self.dqx) + np.square(self.dqy))
        dq = dq[condition]

        return MomentumTransfer._bin_into_q1d(q, dq, i, i_sigma,
                                              self.prefix, bins,
                                              statistic, suffix)

    def bin_annular_into_q1d(self,
                             q_min=0.001,
                             q_max=0.4,
                             bins=100,
                             statistic='mean',
                             suffix="_annular_iq"):
        """
        Wedge calculation and integration

        Calculates: I(Q) and Dq
        The ws_* input parameters are the output workspaces from bin_into_q2d

        Parameters
        ----------
        q_min : float, optional
            , by default
        q_max : float, optional
            , by default
        bins : int or sequence of scalars, optional
            See `scipy.stats.binned_statistic`.
            If `bins` is an int, it defines the number of equal-width bins in
            the given range (10 by default).  If `bins` is a sequence, it
            defines the bin edges, including the rightmost edge, allowing for
            non-uniform bin widths.  Values in `x` that are smaller than lowest
            bin edge areassigned to bin number 0, values beyond the highest bin
            are assigned to ``bins[-1]``.  If the bin edges are specified,
            the number of bins will be, (nx = len(bins)-1).
        statistic : str, optional
            See `scipy.stats.binned_statistic`.
            The statistic to compute, by default 'mean'
            The following statistics are available:
            * 'mean' : compute the mean of values for points within each bin.
                Empty bins will be represented by NaN.
            * 'std' : compute the standard deviation within each bin. This
                is implicitly calculated with ddof=0.
            * 'median' : compute the median of values for points within each
                bin. Empty bins will be represented by NaN.
            * 'count' : compute the count of points within each bin.  This is
                identical to an unweighted histogram.  `values` array is not
                referenced.
            * 'sum' : compute the sum of values for points within each bin.
                This is identical to a weighted histogram.
            * 'min' : compute the minimum of values for points within each bin.
                Empty bins will be represented by NaN.
            * 'max' : compute the maximum of values for point within each bin.
                Empty bins will be represented by NaN.
            * function : a user-defined function which takes a 1D array of
                values, and outputs a single numerical statistic. This function
                will be called on the values in each bin.  Empty bins will be
                represented by function([]), or NaN if this returns an error.
        out_ws_prefix : str, optional
            The prefix of the workspace created in Mantid, by default "ws"

        Returns
        -------
        tuple
            (workspace name, workspace)

        """

        q = np.sqrt(np.square(self.qx) + np.square(self.qy))
        condition = (q >= q_min) & (q <= q_max)

        q = q[condition]
        i = self.i[condition]
        i_sigma = self.i_sigma[condition]

        dq = np.sqrt(np.square(self.dqx) + np.square(self.dqy))
        dq = dq[condition]

        return MomentumTransfer._bin_into_q1d(q, dq, i, i_sigma,
                                              self.prefix, bins,
                                              statistic, suffix)

