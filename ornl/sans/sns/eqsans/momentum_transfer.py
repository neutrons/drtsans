from __future__ import print_function

import numpy as np

from mantid import mtd
from mantid.simpleapi import AddSampleLog, ExtractSpectra, Rebin
from ornl.sans.momentum_transfer import \
    MomentumTransfer as MomentumTransferMain
from ornl.sans.sns.eqsans import geometry
from mantid.kernel import logger
from ornl.sans.samplelogs import SampleLogs


# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide='ignore', invalid='ignore')


class MomentumTransfer(MomentumTransferMain):
    def __init__(self,
                 input_workspace=None,
                 component_name="detector1",
                 out_ws_prefix="ws"):
        super().__init__(input_workspace,
                         component_name=component_name,
                         out_ws_prefix=out_ws_prefix)

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

        if self.component is None:
            self.component = other.component
            self.detector_dims = (self.component.dim_x, self.component.dim_y)

        return self


###############################################################################
# API
###############################################################################


def prepare_momentum_transfer(input_workspace,
                              wavelength_binning=[0.5],
                              sample_aperture=10.0,
                              prefix=None,
                              suffix="_iqxqy_table"):
    """Generates the table workspace necessary for the binning.
    This table contains unbinned Qx Qy.
    It is named `prefix + suffix`, by default:
        `input_workspace.name() + "_iqxqy_table"`

    Parameters
    ----------
    input_workspace : Workspace2D
        The corrected Workspace
    wavelength_binning : list, optional
        use: [min, step, max] or [step]
        By default [0.5], i.e. the step 0.5.
        This is the same as mantid binning.
    sample_aperture : float, optional
        Sample aperture diameter, by default 10.0
    prefix : string, optional
        if None uses `input_workspace.name()`, by default None
    suffix : str, optional
        The suffix for the table workspace, by default "_iqxqy_table"

    Returns
    -------
        Workspace2D
        or
        (Workspace2D, Workspace2D) in case of frame skipping dat aset
    """

    AddSampleLog(Workspace=input_workspace,
                 LogName='sample-aperture-diameter',
                 LogText='{}'.format(sample_aperture),
                 LogType='Number',
                 LogUnit='mm')
    geometry.source_aperture_diameter(input_workspace)

    frames = []
    sl = SampleLogs(input_workspace)
    if bool(sl.is_frame_skipping.value) is True:

        logger.information(
            "This is a frame skipping data set. Outputing 2 WS.")
        frame1_wavelength_min = sl.wavelength_skip_min.value
        frame1_wavelength_max = sl.wavelength_skip_max.value

        frame2_wavelength_min = sl.wavelength_lead_min.value
        frame2_wavelength_max = sl.wavelength_lead_max.value

        frames.append((frame1_wavelength_min, frame1_wavelength_max))
        frames.append((frame2_wavelength_min, frame2_wavelength_max))
    else:
        frame1_wavelength_min = sl.wavelength_min.value
        frame1_wavelength_max = sl.wavelength_max.value
        frames.append((frame1_wavelength_min, frame1_wavelength_max))

    if len(wavelength_binning) > 1 and bool(sl.is_frame_skipping.value) is True:
        logger.error(
            "The WS is frame skipping, use only the step in wavelength_binning")
        return

    for index, (wavelength_min, wavelength_max) in enumerate(frames):

        if len(wavelength_binning) == 1:
            wavelength_binning = [wavelength_min,
                                  wavelength_binning[0], wavelength_max]

        ws_rebin = Rebin(
            InputWorkspace=input_workspace,
            OutputWorkspace="ws_tmp_rebin",
            Params=wavelength_binning)

        bins = ws_rebin.readX(0)
        bin_step = abs(bins[1] - bins[0])

        if prefix is None:
            prefix = input_workspace.name()

        mt_sum = MomentumTransfer(out_ws_prefix=prefix)

        for bin_start in bins[:-1]:
            ws_extracted = ExtractSpectra(InputWorkspace=ws_rebin,
                                          XMin=bin_start,
                                          XMax=bin_start + bin_step)
            wavelength_mean = bin_start + bin_step / 2.0
            AddSampleLog(Workspace=ws_extracted,
                         LogName='wavelength',
                         LogText="{:.2f}".format(wavelength_mean),
                         LogType='Number',
                         LogUnit='Angstrom')
            AddSampleLog(Workspace=ws_extracted,
                         LogName='wavelength-spread',
                         LogText='0.2',
                         LogType='Number',
                         LogUnit='Angstrom')

            mt_extracted = MomentumTransfer(ws_extracted)
            mt_sum += mt_extracted

        if len(frames) == 1:
            return mt_sum.q2d(suffix=suffix)
        else:
            mt_sum.q2d(suffix=suffix+"_frame_{}".format(index+1))

    return mtd[prefix+suffix+"_frame_1"], mtd[prefix+suffix+"_frame_2"]


def iq(input_table_workspace, bins=100, log_binning=False, suffix="_iq"):
    """
    Creates a WS named: input_workspace.name() + "_iq"

    Parameters
    ----------
    input_table_workspace : TableWorkspace

    bins : int or sequence of scalars, optional
        See `scipy.stats.binned_statistic`.
        If `bins` is an int, it defines the number of equal-width bins in
        the given range (10 by default).  If `bins` is a sequence, it
        defines the bin edges, including the rightmost edge, allowing for
        non-uniform bin widths.  Values in `x` that are smaller than lowest
        bin edge areassigned to bin number 0, values beyond the highest bin
        are assigned to ``bins[-1]``.  If the bin edges are specified,
        the number of bins will be, (nx = len(bins)-1).
    log_binning : bool, optional
        if True bins must be an integer, by default False
    """

    mt = MomentumTransfer(input_table_workspace,
                          out_ws_prefix=input_table_workspace.name())

    if log_binning and isinstance(bins, int):
        # TODO: calculate and keep q in MomentumTransfer?
        q = np.sqrt(np.square(mt.qx) + np.square(mt.qy))
        bins = np.logspace(np.log10(np.min(q)), np.log10(np.max(q)), num=bins)

    _, ws = mt.bin_into_q1d(bins=bins, suffix=suffix)
    return ws


def iqxqy(input_table_workspace, bins=100, suffix='_iqxqy'):
    """
    Creates a WS named: input_workspace.name() + "_iqxqy"

    Parameters
    ----------
    input_table_workspace : TableWorkspace

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

    """

    mt = MomentumTransfer(input_table_workspace,
                          out_ws_prefix=input_table_workspace.name())

    # TODO: log binning
    # log_binning : bool, optional
    #     if True bins must be an integer or an array with two integers, e.g.,
    #     100 or [100, 100]
    #     , by default False
    # if log_binning:
    #     num_x = None
    #     num_y = None
    #     if isinstance(bins, int):
    #         num_x = bins
    #         num_y = bins
    #     elif isinstance(bins, list) and len(bins) == 2:
    #         num_x = bins[0]
    #         num_y = bins[1]
    #     if num_x is not None and num_y is not None:
    #         bins_qx = np.logspace(np.log10(np.min(mt.qx)),
    #                               np.log10(np.max(mt.qx)),
    #                               num=num_x)
    #         bins_qy = np.logspace(np.log10(np.min(mt.qy)),
    #                               np.log10(np.max(mt.qy)),
    #                               num=num_y)
    #         bins = [bins_qx, bins_qy]

    _, ws = mt.bin_into_q2d(bins=bins, suffix=suffix)
    return ws
