from enum import Enum
from typing import Any, List, Tuple, Union

import numpy
import numpy as np
from mantid.simpleapi import logger

from drtsans.dataobjects import (
    DataType,
    IQazimuthal,
    IQmod,
    I1DAnnular,
    concatenate,
    getDataType,
    q_azimuthal_to_q_modulo,
)
from drtsans.determine_bins import (
    BinningParams,
    Bins,
    determine_1d_linear_bins,
    determine_1d_log_bins,
)

# To ignore warning:   invalid value encountered in true_divide
np.seterr(divide="ignore", invalid="ignore")

__all__ = [
    "bin_all",
    "bin_intensity_into_q1d",
    "select_i_of_q_by_wedge",
    "bin_annular_into_q1d",
    "bin_intensity_into_q2d",
    "BinningMethod",
    "check_iq_for_binning",
    "determine_1d_linear_bins",
    "determine_1d_log_bins",
    "BinningParams",
]


class BinningMethod(Enum):
    """
    Binning method
    """

    NOWEIGHT = 1  # no-weight binning
    WEIGHTED = 2  # weighted binning


def check_iq_for_binning(i_of_q: Union[IQmod, IQazimuthal]) -> bool:
    """Check I(Q) for binning.

    Binning I(Q) assumes that
    1. there is no NaN or Infinity in intensities
    2. there is no NaN, Infinity or Zero in intensity errors

    Parameters
    ----------
    i_of_q: Union[IQmod, IQazimuthal]

    Returns
    -------
    bool
        True if the input I(Q) for binning meets the assumptions
    """

    error_message = ""

    # Check intensity
    if np.where(np.isnan(i_of_q.intensity))[0].size > 0:
        error_message += "Intensity has {} NaNs: {}\n".format(
            len(np.where(np.isnan(i_of_q.intensity))[0]),
            np.where(np.isnan(i_of_q.intensity))[0],
        )
    if np.where(np.isinf(i_of_q.intensity))[0].size > 0:
        error_message += "Intensity has {} Infinities: {}\n".format(
            len(np.where(np.isnan(i_of_q.intensity))[0]),
            np.where(np.isnan(i_of_q.intensity))[0],
        )

    # Check error
    if np.where(np.isnan(i_of_q.error))[0].size > 0:
        nan_error = np.where(np.isnan(i_of_q.error))[0]
        error_message += f"Intensity error has {len(nan_error)} NaNs: {nan_error}\n"
    if np.where(np.isinf(i_of_q.error))[0].size > 0:
        inf_error = np.where(np.isinf(i_of_q.error))[0]
        error_message += f"Intensity error has {len(inf_error)} Infinities: {inf_error}\n"
    if np.where(np.abs(i_of_q.error) < 1e-20)[0].size > 0:
        zero_error = np.where(np.abs(i_of_q.error) < 1e-20)[0]
        error_message += f"Intensity error has {len(zero_error)} zeros: {zero_error}\n"

    if len(error_message) > 0:
        logger.warning(f"Input I(Q) for binning does not meet assumption: {error_message}")

    return len(error_message) == 0


def valid_wedge(min_angle, max_angle) -> List[Tuple[float, float]]:
    """
    Helper function to validate wedge. It checks that the values  are in the [-90,270) range and
    that the wedge angle is less than 180 degrees

    Parameters
    ----------
    min_angle: float
    max_angle: float

    Returns
    -------
    ~list
        (min_angle, max_angle) tuple. If min_angle < 270 and max_angle > -90
        the function returns two wedges [(min_angle,270.1),(-90.1,max_angle)]
    """
    if min_angle >= 270.0 or min_angle < -90:
        raise ValueError("minimum angle not in the [-90,270) range: {:.1f}".format(min_angle))
    if max_angle >= 270.0 or max_angle < -90:
        raise ValueError("maximum angle not in the [-90,270) range: {:.1f}".format(max_angle))
    if max_angle == min_angle:
        raise ValueError("maximum angle = minimum angle: {:.1f} == {:.1f}".format(min_angle, max_angle))
    if max_angle > min_angle:
        diff = max_angle - min_angle
        if diff < 180.0:
            return [(min_angle, max_angle)]
        raise ValueError(
            "wedge angle is greater than 180 degrees: {:.1f} - {:.1f} = {:.1f} < 180".format(
                max_angle, min_angle, diff
            )
        )
    else:
        # max_angle is smaller than min_angle, meaning that the wedge crosses the 270 -> -90 border
        # Split the wedge into two parts on either side of the border so that all wedges are in the
        # range [-90,270).
        diff = 360.0 - (min_angle - max_angle)
        if diff >= 180:
            raise ValueError(
                f"wedge angle (min {min_angle:.1f} max {max_angle:.1f}) is greater than 180 degrees:"
                f" (360 + {max_angle:.1f}) - {min_angle:.1f} = {diff:.1f} >= 180"
            )
        return [(min_angle, 270.1), (-90.1, max_angle)]


def get_wedges(min_angle: float, max_angle: float, symmetric_wedges=True) -> List[Tuple[float, float]]:
    """
    Helper function to return all wedges defined by the min_angle and max_angle,
    including the wedge offset by 180 degrees

    Parameters
    ----------
    min_angle: float
        lower boundary of the wedge angle in degree
    max_angle: float
        upper boundary of the wedge angle in degree
    symmetric_wedges: bool
        Add the wedge offset by 180 degrees if True

    Returns
    -------
    ~list
        (min_angle, max_angle) tuples.
    """
    if symmetric_wedges:
        wedges = valid_wedge(min_angle, max_angle)
        # create the opposite side and add it
        opp_min = min_angle + 180.0
        opp_max = max_angle + 180.0
        if opp_min >= 270:
            opp_min -= 360.0
        if opp_max >= 270.0:
            opp_max -= 360
        wedges.extend(valid_wedge(opp_min, opp_max))
    elif isinstance(min_angle, (float, int)):
        # also a tuple for a single wedge (min angle, max angle)
        # but not symmetric
        wedges = valid_wedge(min_angle, max_angle)

    else:
        # in this case min_angle and max_angle are actually two wedges
        # that should be summed together
        raise NotImplementedError("use case to have min_angle and max_angle as 2-tuple is disabled")

    return wedges


def bin_all(
    i_qxqy,
    i_modq,
    nxbins,
    nybins,
    n1dbins=None,
    n1dbins_per_decade=None,
    bin1d_type="scalar",
    log_scale=False,
    decade_on_center=False,
    qmin=None,
    qmax=None,
    wedge1_qmin=None,
    wedge1_qmax=None,
    wedge2_qmin=None,
    wedge2_qmax=None,
    qxrange=None,
    qyrange=None,
    annular_angle_bin=1.0,
    wedges: List[Any] = None,
    symmetric_wedges: bool = True,
    error_weighted=False,
    n_wavelength_bin=1,
) -> Tuple[IQazimuthal, List[IQmod]]:
    r"""Do all 1D and 2D binning for a configuration or detector

    Parameters
    ----------
    i_qxqy: ~drtsans.dataobjects.IQazimuthal
        Object containing 2D unbinned data I(Qx, Qy). It will be used for 2D binned data,
        and 1D wedge or annular binned data
    i_modq: ~drtsans.dataobjects.IQmod
        Object containing 1D unbinned data I(\|Q\|). It will be used for scalar binned data
    nxbins: int
        number of bins in the x direction for 2D binning
    nybins: int
        number of bins in the y direction for 2D binning
    n1dbins: int
        number of bins for the 1d binning.
    n1dbins_per_decade: int
        Total number of bins will be this value multiplied by
        number of decades from X min to X max
    bin1d_type: str
        type of binning for 1D data. Possible choices are 'scalar', 'annular', or 'wedge'
    log_scale: bool
        if True, 1D scalar or wedge binning will be logarithmic. Ignored for anything else
    decade_on_center: bool
        Flag to have the min X and max X on bin center; Otherwise, they will be on bin boundary
    qmin: float
        Minimum value of the momentum transfer modulus Q
    qmax: float
        Maximum value of the momentum transfer modulus Q
    wedge1_qmin: float
        Minimum value of the momentum transfer modulus Q for the first wedge when ``bin1d_type = 'wedge'``
    wedge1_qmax: float
        Maximum value of the momentum transfer modulus Q for the first wedge when ``bin1d_type = 'wedge'``
    wedge2_qmin: float
        Minimum value of the momentum transfer modulus Q for the second wedge when ``bin1d_type = 'wedge'``
    wedge2_qmax: float
        Maximum value of the momentum transfer modulus Q for the second wedge when ``bin1d_type = 'wedge'``
    qxrange: ~tuple
        qx min and qx max
    qyrange: ~tuple
        qy min and qy max
    annular_angle_bin: float
        width of annular bin in degrrees. Annular binning is linear
    wedges: list
        list of tuples (angle_min, angle_max) for the wedges. Both numbers have to be in
        the [-90,270) range. It will add the wedge offset by 180 degrees dependent
        on ``symmetric_wedges``
    symmetric_wedges: bool
        It will add the wedge offset by 180 degrees if True
    error_weighted: bool
        if True, the binning is done using the Weighted method
    n_wavelength_bin: None, int
        None: keep original wavelength vector.  int: number of wavelength bins.  1 to sum all

    Returns
    -------
    (~drtsans.dataobjects.IQazimuthal, ~list)
        binned IQazimuthal
        list of binned ~drtsans.dataobjects.IQmod objects. The list has length
        1, unless the 'wedge' mode is selected, when the length is the number of
        original wedges
    """
    # Sanity check
    if n_wavelength_bin is None:
        pass
    elif n_wavelength_bin > 1:
        raise NotImplementedError(f"Case for n_wavelength_bin = {n_wavelength_bin} has not been implemented")

    method = BinningMethod.NOWEIGHT
    if error_weighted:
        method = BinningMethod.WEIGHTED

    # 2D binning
    if qxrange is None:
        # default: data's qx range
        qx_min = np.min(i_qxqy.qx)
        qx_max = np.max(i_qxqy.qx)
    else:
        qx_min, qx_max = qxrange
    binning_x = determine_1d_linear_bins(qx_min, qx_max, nxbins)

    if qyrange is None:
        # default: data's qy range
        qy_min = np.min(i_qxqy.qy)
        qy_max = np.max(i_qxqy.qy)
    else:
        qy_min, qy_max = qyrange
    binning_y = determine_1d_linear_bins(qy_min, qy_max, nybins)

    # bin 2D
    # FIXME - 2D binning does not support weighted binning well.  Force it to no-weight binning
    bin_2d_method = BinningMethod.NOWEIGHT
    binned_q2d = bin_intensity_into_q2d(
        i_qxqy,
        binning_x,
        binning_y,
        method=bin_2d_method,
        wavelength_bins=n_wavelength_bin,
    )

    # 1D binning
    binned_q1d_list = []
    if bin1d_type == "annular":
        # annular binning
        bin_params = BinningParams(0.0, 360.0, int(360.0 / annular_angle_bin))
        kwargs = {"method": method}
        if qmin is not None:
            kwargs["q_min"] = qmin
        if qmax is not None:
            kwargs["q_max"] = qmax
        binned_q1d_list.append(bin_annular_into_q1d(i_qxqy, bin_params, wavelength_bins=n_wavelength_bin, **kwargs))
    else:
        # regular binning including 'scalar' and 'wedge'
        if bin1d_type == "scalar":
            unbinned_1d = [i_modq]
            qmin_list = [qmin] if qmin is not None else [i_modq.mod_q.min()]
            qmax_list = [qmax] if qmax is not None else [i_modq.mod_q.max()]
        elif bin1d_type == "wedge":
            # select Q's by wedge angles
            unbinned_1d = bin_into_wedges(i_qxqy, wedges, symmetric_wedges)
            unbinned_1d_flattened = concatenate(unbinned_1d)
            # there are two wedges with individual (qmin, qmax) limits
            qmin_list = [None, None]
            qmin_list[0] = wedge1_qmin if wedge1_qmin is not None else unbinned_1d_flattened.mod_q.min()
            qmin_list[1] = wedge2_qmin if wedge2_qmin is not None else unbinned_1d_flattened.mod_q.min()
            qmax_list = [None, None]
            qmax_list[0] = wedge1_qmax if wedge1_qmax is not None else unbinned_1d_flattened.mod_q.max()
            qmax_list[1] = wedge2_qmax if wedge2_qmax is not None else unbinned_1d_flattened.mod_q.max()
        else:
            raise ValueError(f"bin1d_type of type {bin1d_type} is not available")

        if log_scale:
            for iub, ub1d in enumerate(unbinned_1d):
                # log bins
                bins_1d = determine_1d_log_bins(
                    qmin_list[iub],
                    qmax_list[iub],
                    n_bins_per_decade=n1dbins_per_decade,
                    n_bins=n1dbins,
                    decade_on_center=decade_on_center,
                )
                # The filter is needed for logarithmic binning so that
                # the qmin and qmax are correctly taken into account
                q_filter = np.where(np.logical_and(ub1d.mod_q >= qmin_list[iub], ub1d.mod_q <= qmax_list[iub]))
                ub1d_filtered = IQmod(
                    ub1d.intensity[q_filter],
                    ub1d.error[q_filter],
                    ub1d.mod_q[q_filter],
                    ub1d.delta_mod_q[q_filter] if ub1d.delta_mod_q is not None else None,
                    ub1d.wavelength[q_filter] if ub1d.wavelength is not None else None,
                )
                binned_q1d = bin_intensity_into_q1d(
                    ub1d_filtered,
                    bins_1d,
                    bin_method=method,
                    wavelength_bins=n_wavelength_bin,
                )
                binned_q1d_list.append(binned_q1d)
        else:
            for iub, ub1d in enumerate(unbinned_1d):
                # linear bins
                bins_1d = determine_1d_linear_bins(qmin_list[iub], qmax_list[iub], n1dbins)
                print(f"Linear binning for {iub}!")
                print(
                    f"Number of NaNs = {np.where(np.isnan(ub1d.intensity))[0].shape}; "
                    f"Number of data points = {ub1d.intensity.shape}"
                )
                binned_q1d_list.append(
                    bin_intensity_into_q1d(
                        ub1d,
                        bins_1d,
                        bin_method=method,
                        wavelength_bins=n_wavelength_bin,
                    )
                )

    return binned_q2d, binned_q1d_list


def bin_into_wedges(i_qxqy, wedges: List[Any], symmetric_wedges: bool) -> List[Any]:
    """

    Parameters
    ----------
    i_qxqy
    wedges
    symmetric_wedges

    Returns
    -------
    list
        list of ~drtsans.dataobjects.IQmod

    """
    unbinned_1d = list()

    # Group is an element of the list
    # Each group is a list of 2-tuples, each is for an individual wedge
    validated_wedge_angles_groups = validate_wedges_groups(wedges, symmetric_wedges)

    # Bin!
    for wedge_angles in validated_wedge_angles_groups:
        # select I(Q) by wedge angles
        wedge_pieces = [select_i_of_q_by_wedge(i_qxqy, min_angle, max_angle) for min_angle, max_angle in wedge_angles]

        # concatenate selected I(Q)
        unbinned_1d.append(q_azimuthal_to_q_modulo(concatenate(wedge_pieces)))

    return unbinned_1d


def validate_wedges_groups(wedges, symmetric_wedges) -> List[List[Tuple[float, float]]]:
    """Validate a list of wedges groups

    Parameters
    ----------
    wedges: ~list
        List of wedges group.  Each wedge group is either a list of 2-tuples or a 2-tuple.  Each 2-tuple is a wedge
    symmetric_wedges: bool
        Flag to include all the symmetric wedges of the given wedges to output

    Returns
    -------
    ~list
        List of wedges group.  Each wedge group is a list of 2-tuples.  Each 2-tupel is a wedge

    """
    validated_wedge_angles_groups = list()

    for wedge in wedges:
        # For a given wedge group, it can be a single wedge or a list of wedges
        if isinstance(wedge, tuple):
            # manual wedge: each wedge group contains 1 and only 1 wedge.
            min_wedge_angle, max_wedge_angle = wedge
            wedge_angles = get_wedges(min_wedge_angle, max_wedge_angle, symmetric_wedges)
        elif isinstance(wedge, list):
            # auto wedge: each wedge group contain 2 wedges
            if len(wedges) != 2 or symmetric_wedges:
                # Note: auto wedge shall not have a pair of wedges sent to this method
                # It is worth to discuss how to work with auto/manual wedge with symmetric/asymmetric combination
                # by unified data structure
                raise NotImplementedError(
                    f"Unsupported scenario for automated wedge: number of wedges {len(wedges)}"
                    f" is not equal to 2.  And/or symmetric wedge option {symmetric_wedges} "
                    f"cannot be True."
                )
            wedge_angles = get_wedges(wedge[0][0], wedge[0][1], symmetric_wedges)
            wedge_angles.extend(get_wedges(wedge[1][0], wedge[1][1], symmetric_wedges))
        else:
            # supported wedge group type
            raise TypeError(f"Wedge group {wedges} of type {type(wedges)} is not supported")

        # add the corrected wedge angles
        validated_wedge_angles_groups.append(wedge_angles)

    return validated_wedge_angles_groups


def bin_intensity_into_q1d(
    i_of_q: IQmod, q_bins: Bins, bin_method: BinningMethod = BinningMethod.NOWEIGHT, wavelength_bins: int = 1
) -> IQmod:
    """Binning I(Q) from scalar Q (1D) with linear binning on Q

    Replace intensity, intensity_error, scalar_q, scalar_dq by IQmod
    Replace bins, q_min=None, q_max=None by BinningParams
    bins: number of bins for linear binning; step per decade for logarithm binning
    q_min : Default to min(scalar_q)
    q_max : Default to max(scalar_q)

    Parameters
    ----------
    i_of_q : ~drtsans.dataobjects.IQmod
        Scalar I(Q) including intensity, intensity_error, scalar_q, scalar_dq in 1d nparray
        including: intensity error mod_q delta_mod_q
    q_bins : Bins
        namedtuple for arbitrary bin edges and bin centers
    bin_method : ~drtsans.BinningMethod
        weighted binning or no-weight binning method
    wavelength_bins: None, int
        number of binned wavelength.  If None, do not bin.  If equal to 1, bin all wavelength together

    Returns
    -------
    drtsans.dataobjects.IQmod
        the one dimensional data as a named tuple
    """
    # Check whether input I(Q) meets assumptions
    check_iq_for_binning(i_of_q)

    # bin I(Q)
    if bin_method == BinningMethod.WEIGHTED:
        # weighed binning
        binned_q = _do_1d_weighted_binning(
            i_of_q.mod_q,
            i_of_q.delta_mod_q,
            i_of_q.intensity,
            i_of_q.error,
            q_bins,
            i_of_q.wavelength,
            wavelength_bins,
        )
    else:
        # no-weight binning
        binned_q = _do_1d_no_weight_binning(
            i_of_q.mod_q,
            i_of_q.delta_mod_q,
            i_of_q.intensity,
            i_of_q.error,
            q_bins,
            i_of_q.wavelength,
            wavelength_bins,
        )

    return binned_q


def select_i_of_q_by_wedge(i_of_q, min_wedge_angle, max_wedge_angle):
    """Select a sub set of I(Q) by 2D wedge

    Parameters
    ----------
    i_of_q : ~drtsans.dataobjects.IQazimuthal
         "intensity": intensity, "error": sigma(I), "qx": qx, "qy": qy, "delta_qx": dqx, "delta_qy", dqy
    min_wedge_angle : float
        minimum value of phi/azimuthal angle for wedge
    max_wedge_angle : float
        maximum value of phi/azimuthal angle for wedge

    Returns
    -------
    drtsans.dataobjects.IQazimuthal
        subset of input I(Qx, Qy) with (Qx, Qy) inside defined wedge

    """
    # Calculate wedge angles for each I(Qx, Qy)
    # calculate azimuthal angles from -180 to 180 degrees
    azimuthal_array = np.arctan2(i_of_q.qy, i_of_q.qx) * 180.0 / np.pi
    # correct azimuthal angles to -90 to 270 degrees
    azimuthal_array[azimuthal_array < -90.0] += 360.0

    # Define the filter (mask/ROI) for pixels falling into preferred wedge
    wedge_indexes = (azimuthal_array > min_wedge_angle) & (azimuthal_array < max_wedge_angle)

    # Create a new IQazimuthal with selected subset
    wedge_i_of_q = IQazimuthal(
        i_of_q.intensity[wedge_indexes],
        i_of_q.error[wedge_indexes],
        i_of_q.qx[wedge_indexes],
        i_of_q.qy[wedge_indexes],
        i_of_q.delta_qx[wedge_indexes],
        i_of_q.delta_qy[wedge_indexes],
    )

    return wedge_i_of_q


def _toQmodAndAzimuthal(
    data: IQazimuthal,
) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """This function returns the values of qmod and azimuthal that are parallel
    to the original data array. It requires that the data is IQazimuthal

    Parameters
    ==========
    data: ~drtsans.dataobjects.IQazimuthal

    Results
    =======
    tuple
        ``(intensity, error, qmod, dqmod, azimuthal)`` as 1d arrays
        with Q in angstrom and azimuthal angle in degrees"""
    if not getDataType(data) == DataType.IQ_AZIMUTHAL:
        raise RuntimeError("Calculating qmod and azimuthal only works for IQazimuthal")

    # reshape the qx and qy if intensity array is 2d
    if len(data.intensity.shape) == 2 and len(data.qx.shape) == 1 and len(data.qy.shape) == 1:
        qx = np.tile(data.qx, (data.qy.shape[0], 1))
        qy = np.tile(data.qy, (data.qx.shape[0], 1)).transpose()
    else:
        qx = data.qx
        qy = data.qy

    # calculate q-scalar
    q = np.sqrt(np.square(qx) + np.square(qy)).ravel()

    # calculate dQ from dQx and dQy
    if data.delta_qx is None or data.delta_qy is None:
        dq = None
    else:
        dq = np.sqrt(np.square(data.delta_qx) + np.square(data.delta_qy)).ravel()

    # azimuthal is expected to be positive so use cyclical nature of trig functions
    azimuthal = np.rad2deg(np.arctan2(qy, qx)).ravel()
    azimuthal[azimuthal < 0.0] += 360.0

    return data.intensity.ravel(), data.error.ravel(), q, dq, azimuthal


def bin_annular_into_q1d(
    i_of_q, phi_bin_params, q_min=0.001, q_max=0.4, method=BinningMethod.NOWEIGHT, wavelength_bins: int | None = 1
):
    """Annular 1D binning

    Calculates: I(azimuthal) and sigma I by assigning pixels to proper azimuthal angle bins
    Input is I(Qx, Qy) and output is :py:obj:`~drtsans.dataobjects.I1DAnnular`.
    The independent axis is the azimuthal angle around the beam center.

    Parameters
    ----------
    i_of_q :  ~drtsans.dataobjects.IQazimuthal
        I(Qx, Qy), sigma I(Qx, Qy), Qx, Qy, dQx and dQy
    phi_bin_params : ~drtsans.BinningParams
        binning parameters on annular angle 'phi'

        phi_min : float
            minimum value of phi/azimuthal angle
        phi_max : float
            maximum value of phi/azimuthal angle
        bins : int or sequence of scalars, optional
            See `scipy.stats.binned_statistic`.
            If `bins` is an int, it defines the number of equal-width bins in
            the given range (10 by default).  If `bins` is a sequence, it
            defines the bin edges, including the rightmost edge, allowing for
            non-uniform bin widths.  Values in `x` that are smaller than lowest
            bin edge areassigned to bin number 0, values beyond the highest bin
            are assigned to ``bins[-1]``.  If the bin edges are specified,
            the number of bins will be, (nx = len(bins)-1).

    q_min : float, optional
        , by default
    q_max : float, optional
        , by default
    method : ~drtsans.BinningMethod
        binning method, no-weight or weighed
    wavelength_bins: None, int
        number of binned wavelength. If None, do not bin. If equal to 1, bin all wavelength together

    Returns
    -------
    I1DAnnular
        Annular-binned I(azimuthal) in 1D
    """
    # Determine azimuthal angle bins (i.e., phi bins)
    phi_bins = determine_1d_linear_bins(phi_bin_params.min, phi_bin_params.max, phi_bin_params.bins)
    if phi_bins.centers.min() < 0.0 or phi_bins.centers.max() > 360.0:
        msg = "must specify range 0<=phi<=360deg found {}<=phi<={}deg".format(phi_bin_params.min, phi_bin_params.max)
        raise ValueError(msg)

    # Check whether input I(Q) meets assumptions
    check_iq_for_binning(i_of_q)

    # convert the data to q and azimuthal angle
    intensity, error, q_array, dq_array, phi_array = _toQmodAndAzimuthal(i_of_q)

    # Filter by q_min and q_max
    allowed_q_index = np.logical_and((q_array > q_min), (q_array < q_max))

    # Construct the wavelength array
    wl_array = None
    if wavelength_bins is None and i_of_q.wavelength is not None:
        wl_array = i_of_q.wavelength[allowed_q_index]

    # select binning method
    # the methods call the independent axis "Q", but are generic to whatever values are passed in
    if method == BinningMethod.NOWEIGHT:
        # no weight binning
        do_1d_binning = _do_1d_no_weight_binning
    elif method == BinningMethod.WEIGHTED:
        # weighted binning
        do_1d_binning = _do_1d_weighted_binning
    else:
        # not supported case
        raise RuntimeError("Binning method {} is not recognized".format(method))

    # apply the selected binning method by either using or skipping the dq_array
    if dq_array is None:
        binned_i_annular = do_1d_binning(
            phi_array[allowed_q_index],
            None,
            intensity[allowed_q_index],
            error[allowed_q_index],
            phi_bins,
            wl_array,
            wavelength_bins,
            annular=True,
        )
    else:
        binned_i_annular = do_1d_binning(
            phi_array[allowed_q_index],
            dq_array[allowed_q_index],
            intensity[allowed_q_index],
            error[allowed_q_index],
            phi_bins,
            wl_array,
            wavelength_bins,
            annular=True,
        )

    return binned_i_annular


def _do_1d_no_weight_binning(
    q_array, dq_array, iq_array, sigmaq_array, q_bins, wl_array, wavelength_bins, annular=False
):
    """Bin I(Q) by given bin edges and do no-weight binning

    This method implements equation 11.34, 11.35 and 11.36 in master document.

    If there is no Q in a certain Qk bin, NaN will be set to both I(Qk) and sigma I(Qk)

    Parameters
    ----------
    q_array: ndarray
        scalar momentum transfer Q in 1D array flattened from 2D detector
    dq_array: ndarray, None
        scalar momentum transfer (Q) resolution in 1D array flattened from 2D detector
    iq_array: ndarray
        I(Q) in 1D array flattened from 2D detector
    sigmaq_array: ndarray
        sigma I(Q) in 1D array flattened from 2D detector
    q_bins: ~drtsans.determine_bins.Bins
        Bin centers and edges
    annular: bool
        If True, returns an I1DAnnular object rather than an IQmod object. Default is False.

    Returns
    -------
    IQmod | I1DAnnular
        IQmod and I1DAnnular are classes for holding 1D binned data.

    """

    def _bin_iq1d(q_bins, q_vec, dq_vec, i_vec, error_vec):
        """Bin I(Q1D), dI(Q1D) and dQ(Q1D) by no weight binning algorithm

        Parameters
        ----------
        q_bins: ~drtsans.determine_bins.Bins
            bin edges and centers
        q_vec: ~numpy.ndarray
            vector of Q1D
        dq_vec: ~numpy.ndarray, None
            vector for Q1D resolution.  could be left as None
        i_vec: ~numpy.ndarray
            1d vector of intensity
        error_vec: ~numpy.ndarray
            1D vector of intensity error

        Returns
        -------
        ~tuple
            binned intensity vector,
            binned intensity error vector,
            binned q resolution vector,
            binned q center vector

        """

        # Count number of Q in 'q_array' in each Q-bin when they are binned (histogram) to 'bin_edges'
        num_pt_vec, _ = np.histogram(q_vec, bins=q_bins.edges)

        # Counts per bin: I_{k, raw} = \sum I(i, j) for each bin
        i_raw_vec, _ = np.histogram(q_vec, bins=q_bins.edges, weights=i_vec)

        # Square of summed uncertainties for each bin
        sigma_sqr_vec, _ = np.histogram(q_vec, bins=q_bins.edges, weights=error_vec**2)

        all_lists = iter([num_pt_vec, i_raw_vec, sigma_sqr_vec])
        assert all(len(lst) == len(num_pt_vec) for lst in all_lists)

        # Final I(Q):     I_k = \frac{I_{k, raw}}{N_k}
        i_final_vec = i_raw_vec / num_pt_vec

        # Final sigma(Q): sigmaI_k = \frac{sigmaI_{k, raw}}{N_k}
        sigma_final_vec = np.sqrt(sigma_sqr_vec) / num_pt_vec
        # Replace 0 error with 1, as uncertainty is not zero in the region of interest
        sigma_final_vec[sigma_final_vec == 0] = 1.0

        # Calculate Q resolution of binned
        if dq_vec is None:
            bin_dq_vec = None
        else:
            binned_vec, _ = np.histogram(q_vec, bins=q_bins.edges, weights=dq_vec)
            bin_dq_vec = binned_vec / num_pt_vec
            assert len(bin_dq_vec) == len(i_final_vec)

        return i_final_vec, sigma_final_vec, bin_dq_vec

    # check input
    assert q_bins.centers.shape[0] + 1 == q_bins.edges.shape[0]

    if wavelength_bins == 1 or wl_array is None:
        # bin I(Q, wl) regardless of wl value
        binned_i_vec, binned_sigma_vec, binned_dq_vec = _bin_iq1d(q_bins, q_array, dq_array, iq_array, sigmaq_array)
        binned_q_vec = q_bins.centers
        # construct output without wavelength vector
        binned_wl_vec = None

    elif wavelength_bins is None:
        # bin I(Q) with same value of wavelength
        unique_wl_vec = np.unique(wl_array)
        unique_wl_vec.sort()

        # construct a 2D array for filtering
        if dq_array is None:
            wl_matrix = np.array([wl_array, q_array, iq_array, sigmaq_array])
        else:
            wl_matrix = np.array([wl_array, q_array, iq_array, sigmaq_array, dq_array])
        wl_matrix = wl_matrix.transpose()

        # define output
        binned_q_vec = binned_dq_vec = binned_i_vec = binned_sigma_vec = binned_wl_vec = np.ndarray(
            shape=(0,), dtype=float
        )

        for wl_i in unique_wl_vec:
            # filter
            filtered_matrix = wl_matrix[wl_matrix[:, 0] == wl_i]

            # special work with q resolution
            if dq_array is None:
                dq_array_i = None
            else:
                dq_array_i = filtered_matrix[:, 4]

            # bin by Q1D
            i_final_array, sigma_final_array, bin_q_resolution = _bin_iq1d(
                q_bins,
                filtered_matrix[:, 1],
                dq_array_i,
                filtered_matrix[:, 2],
                filtered_matrix[:, 3],
            )

            # build up the final output
            binned_q_vec = np.concatenate((binned_q_vec, q_bins.centers))
            binned_i_vec = np.concatenate((binned_i_vec, i_final_array))
            binned_sigma_vec = np.concatenate((binned_sigma_vec, sigma_final_array))
            if dq_array is not None:
                binned_dq_vec = np.concatenate((binned_dq_vec, bin_q_resolution))
            binned_wl_vec = np.concatenate((binned_wl_vec, np.zeros_like(i_final_array) + wl_i))
        # END-FOR (wl_i)

        if dq_array is None:
            binned_dq_vec = None
    else:
        raise RuntimeError(f"Number of wavlength bins = {wavelength_bins} is not supported")

    # Get the final result by constructing an IQmod or I1DAnnular object defined in ~drtsans.dataobjects.
    # IQmod and I1DAnnular are classes for holding 1D binned data.
    if annular:
        binned_i1d = I1DAnnular(
            intensity=binned_i_vec,
            error=binned_sigma_vec,
            phi=binned_q_vec,
            wavelength=binned_wl_vec,
        )
    else:
        binned_i1d = IQmod(
            intensity=binned_i_vec,
            error=binned_sigma_vec,
            mod_q=binned_q_vec,
            delta_mod_q=binned_dq_vec,
            wavelength=binned_wl_vec,
        )

    return binned_i1d


def _do_1d_weighted_binning(
    q_array, dq_array, iq_array, sigma_iq_array, q_bins, wl_array, wavelength_bins, annular=False
):
    """Bin I(Q) by given bin edges and do weighted binning

    This method implements equation 11.22, 11.23 and 11.24 in master document for 1-dimensional Q

    If there is no Q in a certain Qk bin, NaN will be set to both I(Qk) and sigma I(Qk)

    General description of algorithm:

    Equation 11.26
    I(Q') = sum_{Q, lambda}^{K} (I(Q, lambda) / sigma(Q, lambda)^2) /
            sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)

    Equation 11.27
    sigmaI(Q') = sqrt(sum_{Q, lambda}^{K} (sigma(Q, lambda / sigma(Q, lambda)^2)^2) /
                 sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)
               = sqrt(sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)) /
                 sum_{Q, lambda}^{K}(1/sigma(Q, lambda)^2)
               = 1 / sqrt(sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2))

    Equation 11.28
    sigmaQ(Q') = sum_{Q, lambda}^{K}(sigmaQ(Q, lambda)/sigma^2(Q, lambda)^2) /
                 sum_{Q, lambda}^{K}(1/sigma(Q, lambda)^2)

    Parameters
    ----------
    q_array: ndarray
        scalar momentum transfer Q in 1D array flattened from 2D detector
    dq_array: ndarray, None
        scalar momentum transfer (Q) resolution in 1D array flattened from 2D detector
    iq_array: ndarray
        I(Q) in 1D array flattened from 2D detector
    sigma_iq_array: ndarray
        sigma I(Q) in 1D array flattened from 2D detector
    q_bins: ~drtsans.determine_bins.Bins
        Bin centers and edges
    annular: bool
        If True, returns an I1DAnnular object rather than an IQmod object. Default is False.

    Returns
    -------
    IQmod | I1DAnnular
        IQmod and I1DAnnular are classes for holding 1D binned data.

    """

    def _bin_q1d_weighted(q_bins, mod_q_array, delta_q_array, intensity_array, error_array):
        """Do 1D weighed binning"""
        # bin I(Q, wl) regardless of wl value
        # Calculate 1/sigma^2 for multiple uses
        invert_sigma2_array = 1.0 / (error_array**2)

        # Histogram on 1/sigma^2, i.e., nominator part in Equation 11.22, 11.23 and 11.24
        # sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)
        w_array, _ = np.histogram(mod_q_array, bins=q_bins.edges, weights=invert_sigma2_array)

        # Calculate Equation 11.26: I(Q)
        #  I(Q') = sum_{Q, lambda}^{K} (I(Q, lambda) / sigma(Q, lambda)^2) /
        #              sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)
        # denominator in Equation 11.22: sum_{Q, lambda}^{K} (I(Q, lambda) / sigma(Q, lambda)^2)
        i_raw_array, _ = np.histogram(mod_q_array, bins=q_bins.edges, weights=intensity_array * invert_sigma2_array)

        assert len(i_raw_array) == len(w_array)

        # numerator divided by denominator (11.26)
        binned_intensity_array = i_raw_array / w_array

        # Calculate equation 11.27: sigmaI(Q)
        # sigmaI(Q') = sqrt(sum_{Q, lambda}^{K} (sigma(Q, lambda / sigma(Q, lambda)^2)^2) /
        #              sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)
        #            = sqrt(sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2)) /
        #              sum_{Q, lambda}^{K}(1/sigma(Q, lambda)^2)
        #             = 1 / sqrt(sum_{Q, lambda}^{K} (1 / sigma(Q, lambda)^2))
        # Thus histogrammed sigmaI can be obtained from histogrammed invert_sigma2_array directly
        binned_error_array = 1 / np.sqrt(w_array)

        # Calculate equation 11.28:  sigmaQ (i.e., Q resolution)
        # sigmaQ(Q') = sum_{Q, lambda}^{K}(sigmaQ(Q, lambda)/sigma^2(Q, lambda)^2) /
        #              sum_{Q, lambda}^{K}(1/sigma(Q, lambda)^2)
        # denominator in Equation 11.28: sum_{Q, lambda}^{K}(1/sigma(Q, lambda)^2)
        if delta_q_array is None:
            binned_dq_array = None
        else:
            binned_dq, _ = np.histogram(mod_q_array, bins=q_bins.edges, weights=delta_q_array * invert_sigma2_array)
            assert len(binned_dq) == len(binned_intensity_array)
            # numerator divided by denominator (11.28)
            binned_dq_array = binned_dq / w_array

        return binned_intensity_array, binned_error_array, binned_dq_array

    # Check input
    assert q_bins.centers.shape[0] + 1 == q_bins.edges.shape[0]

    if wl_array is None or wavelength_bins == 1:
        # bin I(Q, wl) regardless of wl value
        binned_i_vec, binned_sigma_vec, binned_dq_vec = _bin_q1d_weighted(
            q_bins, q_array, dq_array, iq_array, sigma_iq_array
        )
        binned_q_vec = q_bins.centers
        # construct output without wavelength vector
        binned_wl_vec = None

    elif wavelength_bins is None:
        # bin I(Q, wl) per wavelength
        unique_wl_vec = np.unique(wl_array)
        unique_wl_vec.sort()

        # construct a 2D array for filtering
        if dq_array is None:
            wl_matrix = np.array([wl_array, q_array, iq_array, sigma_iq_array])
        else:
            wl_matrix = np.array([wl_array, q_array, iq_array, sigma_iq_array, dq_array])
        wl_matrix = wl_matrix.transpose()

        # define output
        binned_q_vec = binned_dq_vec = binned_i_vec = binned_sigma_vec = binned_wl_vec = np.ndarray(
            shape=(0,), dtype=float
        )

        for wl_i in unique_wl_vec:
            # filter
            filtered_matrix = wl_matrix[wl_matrix[:, 0] == wl_i]

            # special work with q resolution
            if dq_array is None:
                dq_array_i = None
            else:
                dq_array_i = filtered_matrix[:, 4]

            # bin by Q1D
            i_final_array, sigma_final_array, bin_q_resolution = _bin_q1d_weighted(
                q_bins=q_bins,
                mod_q_array=filtered_matrix[:, 1],
                delta_q_array=dq_array_i,
                intensity_array=filtered_matrix[:, 2],
                error_array=filtered_matrix[:, 3],
            )

            # build up the final output
            binned_q_vec = np.concatenate((binned_q_vec, q_bins.centers))
            binned_i_vec = np.concatenate((binned_i_vec, i_final_array))
            binned_sigma_vec = np.concatenate((binned_sigma_vec, sigma_final_array))
            if dq_array is not None:
                binned_dq_vec = np.concatenate((binned_dq_vec, bin_q_resolution))
            binned_wl_vec = np.concatenate((binned_wl_vec, np.zeros_like(i_final_array) + wl_i))
        # END-FOR (wl_i)

        if dq_array is None:
            binned_dq_vec = None

    else:
        raise RuntimeError(f"Binning with wavelength bins = {wavelength_bins} is not supported")

    # Get the final result by constructing an IQmod or I1DAnnular object defined in ~drtsans.dataobjects.
    # IQmod and I1DAnnular are classes for holding 1D binned data.
    if annular:
        binned_i1d = I1DAnnular(
            intensity=binned_i_vec,
            error=binned_sigma_vec,
            phi=binned_q_vec,
            wavelength=binned_wl_vec,
        )
    else:
        binned_i1d = IQmod(
            intensity=binned_i_vec,
            error=binned_sigma_vec,
            mod_q=binned_q_vec,
            delta_mod_q=binned_dq_vec,
            wavelength=binned_wl_vec,
        )

    return binned_i1d


def bin_intensity_into_q2d(i_of_q, qx_bins, qy_bins, method=BinningMethod.NOWEIGHT, wavelength_bins=1):
    """Bin I(Qx, Qy) into to new (Qx, Qy) bins

    Note 1: for binning parameters:
    - 'min': float or None.  If None, set to default as min(Qx) (or Qy)
    - 'max': float or None.  If None, set to default as max(Qx) (or Qy)
    - 'bins': integer as number of bins

    Note 2: output Intensity, error, dqx an dqy are in following order
    -    qx = [[qx0, qx1, ...], [qx0, qx1, ...], ...]
    -    qy = [[qy0, qy0, ...], [qy1, qy1, ...], ...]

    Parameters
    ----------
    i_of_q: ~drtsans.dataobjects.IQazimuthal
        class IQazimuthal(namedtuple('IQazimuthal', 'intensity error qx qy delta_qx delta_qy wavelength'))
    qx_bins : Bins
        namedtuple for arbitrary bin edges and bin centers for Qx
    qy_bins : Bins
        namedtuple for arbitrary bin edges and bin centers for Qy
    method: ~drtsans.BinningMethod
        Weighted binning or no weight binning
    wavelength_bins: None, int
        number of binned wavelength.  If None, do not bin.  If equal to 1, bin all wavelength together

    Returns
    -------
    ~drtsans.dataobjects.IQazimuthal
        binned IQazimuthal (important: must read Note 2)
    """
    # Check whether input I(Q) meets assumptions
    check_iq_for_binning(i_of_q)

    # Check whether it needs to bin wavelength
    if wavelength_bins == 1 or i_of_q.wavelength is None:
        # no need to do 2D binning by filtering wavelength
        filter_wavelength = False
    else:
        # bin I(qx, qy, wl) by each wave length
        filter_wavelength = True

    if method == BinningMethod.NOWEIGHT:
        # Calculate no-weight binning
        binned_arrays = _do_2d_no_weight_binning(
            i_of_q.qx,
            i_of_q.delta_qx,
            i_of_q.qy,
            i_of_q.delta_qy,
            i_of_q.wavelength,
            i_of_q.intensity,
            i_of_q.error,
            qx_bins.edges,
            qy_bins.edges,
            debug_filter_wl=filter_wavelength,
            sum_all_wavelengths=not filter_wavelength,
        )
    else:
        # Calculate weighed binning
        # FIXME - _do_2d_weighted_binning API shall be changed as _do_2d_no_weight_binning() for sum_all_wavelength
        binned_arrays = _do_2d_weighted_binning(
            i_of_q.qx,
            i_of_q.delta_qx,
            i_of_q.qy,
            i_of_q.delta_qy,
            i_of_q.wavelength,
            i_of_q.intensity,
            i_of_q.error,
            qx_bins.edges,
            qy_bins.edges,
        )
    # END-IF-ELSE

    # construct return
    binned_intensities, binned_sigmas, binned_dqx, binned_dqy, binned_wl = binned_arrays
    # create Qx and Qy meshgrid explicitly
    # this must agree with the return from histogram2D, which is as
    # qx = [[qx0, qx0, ...],
    #       [qx1, qx1, ...],
    #       ...]
    # qy = [[qy0, qy1, ...],
    #       [qy0, qy1, ...],
    #       ...]
    # Thus indexing='ij' is used
    qx_matrix, qy_matrix = np.meshgrid(qx_bins.centers, qy_bins.centers, indexing="ij")
    binned_qx_array = qx_matrix
    binned_qy_array = qy_matrix
    unique_wl_vec = np.unique(i_of_q.wavelength)
    unique_wl_vec.sort()
    if i_of_q.wavelength is not None and filter_wavelength:
        for wl_i in unique_wl_vec[1:]:
            binned_qx_array = np.concatenate((binned_qx_array, qx_matrix), axis=1)
            binned_qy_array = np.concatenate((binned_qy_array, qy_matrix), axis=1)
    return IQazimuthal(
        intensity=binned_intensities,
        error=binned_sigmas,
        qx=binned_qx_array,
        delta_qx=binned_dqx,
        qy=binned_qy_array,
        delta_qy=binned_dqy,
        wavelength=binned_wl,
    )


def _bin_iq2d(qx_bin_edges, qy_bin_edges, qx_vec, qy_vec, dqx_vec, dqy_vec, i_vec, error_vec):
    """Bin I(Q2D), dI(Q2D) and dQ(Q2D) by no weight binning algorithm

    Parameters
    ----------
    qx_bin_edges: ~numpy.ndarray
        bin edges
    qy_bin_edges: ~numpy.ndarray
        bin edges
    qx_vec: ~numpy.ndarray
        array of Q2D
    qy_vec: ~numpy.ndarray
        array of Q2D
    dqx_vec: ~numpy.ndarray or None
        array for Q2D resolution. May be None
    dqy_vec: ~numpy.ndarray or None
        array for Q2D resolution. May be None
    i_vec: ~numpy.ndarray
        2D array of intensity
    error_vec: ~numpy.ndarray
        2D array of intensity error

    Returns
    -------
    ~tuple
        binned intensity vector, binned intensity error vector,
        binned qx resolution vector, binned qy resolution vector
    """

    # Number of I(q) in each target Q bin
    num_pt_array, *_ = np.histogram2d(qx_vec, qy_vec, bins=(qx_bin_edges, qy_bin_edges))

    # Counts per bin: I_{k, raw} = \sum I(i, j) for each bin
    i_raw_array, *_ = np.histogram2d(qx_vec, qy_vec, bins=(qx_bin_edges, qy_bin_edges), weights=i_vec)

    # Square of summed uncertainties for each bin
    sigma_sqr_array, *_ = np.histogram2d(qx_vec, qy_vec, bins=(qx_bin_edges, qy_bin_edges), weights=error_vec**2)

    # Q resolution: simple average
    dqx_raw_array, *_ = np.histogram2d(qx_vec, qy_vec, bins=(qx_bin_edges, qy_bin_edges), weights=dqx_vec)
    dqy_raw_array, *_ = np.histogram2d(qx_vec, qy_vec, bins=(qx_bin_edges, qy_bin_edges), weights=dqy_vec)

    # Final I(Q): I_{k, final} = \frac{I_{k, raw}}{Nk}
    #       sigma = 1/sqrt(w_k)
    i_final_array = i_raw_array / num_pt_array
    sigma_final_array = np.sqrt(sigma_sqr_array) / num_pt_array
    dqx_final_array = dqx_raw_array / num_pt_array
    dqy_final_array = dqy_raw_array / num_pt_array

    return i_final_array, sigma_final_array, dqx_final_array, dqy_final_array


def _do_2d_no_weight_binning(
    qx_array,
    dqx_array,
    qy_array,
    dqy_array,
    wl_array,
    iq_array,
    sigma_iq_array,
    qx_bin_edges,
    qy_bin_edges,
    sum_all_wavelengths: bool = True,
    debug_filter_wl: bool = False,
):
    """Perform 2D no-weight binning on I(Qx, Qy)

    General description of the algorithm:

      I_{i, j} = sum^{(i, j)}_k I_{k} / N_{i, j}
      sigma I_{i, j} = sqrt(sum^{(i, j)}_k sigma I_k^2) / N_{i, j}

    Parameters
    ----------
    qx_array: ndarray
        Qx array
    dqx_array: ndarray
        Qx resolution
    qy_array : ndarray
        Qy array
    dqy_array: ndarray
        Qy resolution
    wl_array: ndarray or None
        wavelengths
    iq_array: ndarray
        intensities
    sigma_iq_array: ndarray
        intensities error
    qx_bin: ~drtsans.determine_bins.Bins
        Bin centers and edges
    qy_bin:
        Bin centers and edges
    sum_all_wavelengths: bool
        Flag to bin I(qx, qy, wavelength) by qx and qy only.  Wavelength term will then be thrown away

    Returns
    -------
    ndarray, ndarray, ndarray, ndarray, ndarray
        intensities (n x m x o), sigma intensities (n x m x o), Qx resolution (n x m x o), Qy resolution (n x m x o),
        Wavelengths (o)
    """
    if wl_array is None or sum_all_wavelengths:
        # bin only by (qx, qy).  all I(qx, qy, wavelength) with binned regardless of wavelength value
        # output will be I(Qx, Qy)
        (
            binned_iq_array,
            binned_sigma_iq_array,
            binned_dqx_array,
            binned_dqy_array,
        ) = _bin_iq2d(
            qx_bin_edges,
            qy_bin_edges,
            qx_array,
            qy_array,
            dqx_array,
            dqy_array,
            iq_array,
            sigma_iq_array,
        )
        binned_wl_array = None
    else:
        # separate I(qx, qy, wavelength) by wavelength value and bin (qx, qy)
        # output will be I(Qx, Qy, wavelength)
        if debug_filter_wl is False and len(wl_array) > 1:
            raise RuntimeError("It is not supposed to do binning with wavelength term kept.")

        unique_wl_vec = np.unique(wl_array)
        unique_wl_vec.sort()

        # construct a 2D array for filtering
        if dqx_array is None:
            wl_matrix = np.array([wl_array, qx_array, qy_array, iq_array, sigma_iq_array])
        else:
            wl_matrix = np.array(
                [
                    wl_array,
                    qx_array,
                    qy_array,
                    iq_array,
                    sigma_iq_array,
                    dqx_array,
                    dqy_array,
                ]
            )
        wl_matrix = wl_matrix.transpose()

        binned_iq_array = np.ndarray(shape=(0,), dtype=float)
        binned_sigma_iq_array = binned_wl_array = np.ndarray(shape=(0,), dtype=float)

        # Initialize binned dqx and dqy arrays
        if dqx_array is not None:
            binned_dqx_array = np.ndarray(shape=(0,), dtype=float)
            binned_dqy_array = np.ndarray(shape=(0,), dtype=float)
        else:
            binned_dqx_array = binned_dqy_array = None

        for wl_i in unique_wl_vec:
            filtered_matrix = wl_matrix[wl_matrix[:, 0] == wl_i]

            # special work with q resolution
            if dqx_array is None:
                dqx_array_i = None
                dqy_array_i = None
            else:
                dqx_array_i = filtered_matrix[:, 5]
                dqy_array_i = filtered_matrix[:, 6]

            # bin by Q2D
            (
                i_final_array,
                sigma_final_array,
                dqx_final_array,
                dqy_final_array,
            ) = _bin_iq2d(
                qx_bin_edges,
                qy_bin_edges,
                filtered_matrix[:, 1],
                filtered_matrix[:, 2],
                dqx_array_i,
                dqy_array_i,
                filtered_matrix[:, 3],
                filtered_matrix[:, 4],
            )
            # build up the final output
            binned_iq_array = (
                np.concatenate((binned_iq_array, i_final_array), axis=1) if binned_iq_array.size else i_final_array
            )
            binned_sigma_iq_array = (
                np.concatenate((binned_sigma_iq_array, sigma_final_array), axis=1)
                if binned_sigma_iq_array.size
                else sigma_final_array
            )
            if dqx_array is not None:
                binned_dqx_array = (
                    np.concatenate((binned_dqx_array, dqx_final_array), axis=1)
                    if binned_dqx_array.size
                    else dqx_final_array
                )
                binned_dqy_array = (
                    np.concatenate((binned_dqy_array, dqy_final_array), axis=1)
                    if binned_dqy_array.size
                    else dqy_final_array
                )
            binned_wl_array = (
                np.concatenate((binned_wl_array, np.zeros_like(i_final_array) + wl_i), axis=1)
                if binned_wl_array.size
                else np.zeros_like(i_final_array) + wl_i
            )
        # END-FOR (wl_i)

        if dqx_array is None:
            binned_dqx_array = None
            binned_dqy_array = None

    return (
        binned_iq_array,
        binned_sigma_iq_array,
        binned_dqx_array,
        binned_dqy_array,
        binned_wl_array,
    )


def _do_2d_weighted_binning(
    qx_array,
    dqx_array,
    qy_array,
    dqy_array,
    wl_array,
    iq_array,
    sigma_iq_array,
    x_bin_edges,
    y_bin_edges,
):
    """Perform 2D weighted binning

    General description of algorithm:

    I(x', y')      = sum_{x, y, lambda}^{K} (I(x, y, lambda) / sigma(x, y, lambda)^2) /
                     sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
    sigmaI(x', y') = sqrt(sum_{x, y, lambda}^{K} (sigma(x, y, lambda / sigma(x, y, lambda)^2)^2) /
                     sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
                   = sqrt(sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)) /
                     sum_{x, y, lambda}^{K}(1/sigma(x, y, lambda)^2)
                   = 1 / sqrt(sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2))
    sigmaQ(x', y') = sum_{x, y, lambda}^{K}(sigmaQ(x, y, lambda)/sigma^2(x, y, lambda)^2) /
                     sum_{x, y, lambda}^{K}(1/sigma(x, y, lambda)^2)

    where K is the set of (x, y, sigma) such that (x, y, sigma) in the same Q_bin

    Parameters
    ----------
    qx_array : ndarray
        qx
    dqx_array : ndarray
        Qx resolution
    qy_array: ndarray
        qy
    dqy_array: ndarray
        Qy resolution
    wl_array : ndarray
        wavelengths
    iq_array : ndarray
        intensities
    sigma_iq_array : ndarray
        intensity errors
    x_bin_edges : ndarray
        X bin edges
    y_bin_edges
        Y bin edges

    Returns
    -------
    ndarray, ndarray, ndarray, ndarray
        binned intensities (n x m), binned sigmas (n x m), binned Qx resolution (n x m), binned Qy resolution (n x m),
        binned wavelength (n x m)
    """
    unique_wl_vec = np.unique(wl_array)
    unique_wl_vec.sort()
    if wl_array is None or len(unique_wl_vec) == 1:
        # Calculate 1/sigma^2 for multiple uses
        invert_sigma2_array = 1.0 / (sigma_iq_array**2)  # 1D

        # Histogram on 1/sigma^2, i.e., nominator part in Equation 11.22, 11.23 and 11.24
        # sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
        w_2d_array, *_ = np.histogram2d(
            qx_array,
            qy_array,
            bins=(x_bin_edges, y_bin_edges),
            weights=invert_sigma2_array,
        )  # 2D

        # Calculate Equation 11.22: I(Qx, Qy)
        # I(x', y') = sum_{x, y, lambda}^{K} (I(x, y, lambda) / sigma(x, y, lambda)^2) /
        #             sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
        # denominator in Equation 11.22: sum_{x, y, lambda}^{K} (I(x, y, lambda) / sigma(x, y, lambda)^2)
        i_raw_2d_array, *_ = np.histogram2d(
            qx_array,
            qy_array,
            bins=(x_bin_edges, y_bin_edges),
            weights=iq_array * invert_sigma2_array,
        )  # 2D
        # denominator divided by nominator (11.22)
        i_final_array = i_raw_2d_array / w_2d_array

        # Calculate equation 11.23: sigmaI(Q)
        # sigmaI(x', y') = sqrt(sum_{x, y, lambda}^{K} (sigma(x, y, lambda / sigma(x, y, lambda)^2)^2) /
        #                  sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)
        #                = sqrt(sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2)) /
        #                sum_{x, y, lambda}^{K}(1/sigma(x, y, lambda)^2)
        #                = 1 / sqrt(sum_{x, y, lambda}^{K} (1 / sigma(x, y, lambda)^2))
        # Thus histogrammed sigmaI can be obtained from histogrammed invert_sigma2_array directly
        sigma_final_array = 1.0 / np.sqrt(w_2d_array)

        # Calculate equation 11.24:  sigmaQx and sigmaQy (i.e., Q resolution)
        # sigmaQ(x', y') = sum_{x, y, lambda}^{K}(sigmaQ(x, y, lambda)/sigma^2(x, y, lambda)^2) /
        #                  sum_{x, y, lambda}^{K}(1/sigma(x, y, lambda)^2)
        # denominator in Equation 11.24: sum_{x, y, lambda}^{K}(sigmaQ(x, y, lambda)/sigma^2(x, y, lambda)^2)
        dqx_raw_array, *_ = np.histogram2d(
            qx_array,
            qy_array,
            bins=(x_bin_edges, y_bin_edges),
            weights=dqx_array * invert_sigma2_array,
        )  # 2D
        dqy_raw_array, *_ = np.histogram2d(
            qx_array,
            qy_array,
            bins=(x_bin_edges, y_bin_edges),
            weights=dqy_array * invert_sigma2_array,
        )  # 2D
        # denominator divided by nominator (11.24)
        dqx_final_array = dqx_raw_array / w_2d_array  # dQx
        dqy_final_array = dqy_raw_array / w_2d_array  # dQy

        if wl_array is None:
            wl_final_array = None
        else:
            wl_final_array = np.full_like(i_final_array, unique_wl_vec[0])
    else:
        raise NotImplementedError("2D binning with multiple wavelengths is not supported")

    return (
        i_final_array,
        sigma_final_array,
        dqx_final_array,
        dqy_final_array,
        wl_final_array,
    )


def explore_binning_issue(ub_index, n_wavelength_bin, ub1d: IQmod, bins_1d):
    # TODO FIXME - Remove after binning issue is resolved completely

    # DEBUG BINNING
    print(f"[PROOF] [{ub_index}]  Wavelength bins = {n_wavelength_bin}")
    if n_wavelength_bin is None:
        wl_vec = np.unique(ub1d.wavelength)
        print(f"number of wavelength: {len(wl_vec)}: {wl_vec}")
    else:
        wl_vec = None

    for i_bin in range(5):
        # print(f'{i_bin}-bin:  boundary: {bins_1d.edges[i_bin]}, {bins_1d.edges[i_bin + 1]}')
        bin_qmin = bins_1d.edges[i_bin]
        bin_qmax = bins_1d.edges[i_bin + 1]
        # filter the I(Q) in boundary
        # >= q_min
        in_range_i_arrays = ub1d.intensity[ub1d.mod_q >= bin_qmin]
        in_range_q_arrays = ub1d.mod_q[ub1d.mod_q >= bin_qmin]
        # < qmax
        in_range_i_arrays = in_range_i_arrays[in_range_q_arrays < bin_qmax]
        in_range_q_arrays = in_range_q_arrays[in_range_q_arrays < bin_qmax]
        sum_intensity = in_range_i_arrays.sum()
        print(
            f"{i_bin}-bin:  boundary: {bins_1d.edges[i_bin]}, {bins_1d.edges[i_bin + 1]}:"
            f"num points = {len(in_range_q_arrays)}, "
            f"sum = {sum_intensity},"
            f"average = {sum_intensity / len(in_range_q_arrays)}"
        )

        if wl_vec is not None:
            # wavelength in details
            in_range_wl_arrays = ub1d.wavelength[ub1d.mod_q >= bin_qmin]
            in_range_q_arrays = ub1d.mod_q[ub1d.mod_q >= bin_qmin]
            # < qmax
            in_range_wl_arrays = in_range_wl_arrays[in_range_q_arrays < bin_qmax]

            sum_i_per_wl_vec = list()
            num_pt_per_wl_vec = list()
            for wl in wl_vec:
                selected_i_array = in_range_i_arrays[np.abs(in_range_wl_arrays - wl) < 0.001]
                sum_i_per_wl_vec.append(np.sum(selected_i_array))
                num_pt_per_wl_vec.append(len(selected_i_array))
                print(f"wl = {wl}: sum = {np.sum(selected_i_array)}, num I(Q) = {len(selected_i_array)}")
            sum_i_per_wl_vec = np.array(sum_i_per_wl_vec)
            num_pt_per_wl_vec = np.array(num_pt_per_wl_vec)
            sum_i = sum_i_per_wl_vec.sum()

            # binned value can be a little more complicated
            valid_sum_i_vec = sum_i_per_wl_vec[num_pt_per_wl_vec > 0]
            valid_num_pt_vec = num_pt_per_wl_vec[num_pt_per_wl_vec > 0]
            num_valid_ws = len(valid_num_pt_vec)
            bin_int = np.sum(valid_sum_i_vec / valid_num_pt_vec) / num_valid_ws

            print(f"Number of I(Q) = {num_pt_per_wl_vec.sum()}, Total I(Q) = {sum_i}, Binned I = {bin_int}")
