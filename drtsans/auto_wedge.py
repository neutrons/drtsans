import numpy as np
from drtsans.dataobjects import IQmod
from drtsans.determine_bins import determine_1d_linear_bins
from drtsans.iq import BinningMethod, BinningParams, bin_annular_into_q1d
from drtsans.settings import unique_workspace_dundername
# https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
from mantid.simpleapi import DeleteWorkspace, logger
# https://docs.mantidproject.org/nightly/algorithms/Fit-v1.html
from mantid.simpleapi import Fit
import h5py

__all__ = ['getWedgeSelection']

# factor to convert from a sigma to full width half max
# https://docs.mantidproject.org/nightly/fitting/fitfunctions/Gaussian.html
SIGMA_TO_FWHM = 2. * np.sqrt(2. * np.log(2.))


def _debug_output(iq2d, rings, azimuthal_delta):
    """write annular binned I(Qx, Qy) and intensity rings to hdf5 for further review

    Parameters
    ----------
    iq2d:  ~drtsans.dataobjects.Azimuthal
    rings: ~list
        List of I(Qx, Qy) in rings

    Returns
    -------

    """
    # annular binning on the full range
    azimuthal_offset = 0.5 * azimuthal_delta
    azimuthal_binning = BinningParams(0. - azimuthal_offset, 360. - azimuthal_offset,
                                      bins=int(360. / azimuthal_delta))
    q1d = np.sqrt(iq2d.qx**2 + iq2d.qy**2)
    logger.notice(f'[DEBUG] I(Qx, Qy) Q range: {q1d.min()}, {q1d.max()}')
    full_range_annular = bin_annular_into_q1d(iq2d, azimuthal_binning, q1d.min(), q1d.max(),
                                              BinningMethod.NOWEIGHT)

    # open
    debug_h5 = h5py.File('annular_i_q.h5', 'w')

    # write full range
    group = debug_h5.create_group('full range')
    group.create_dataset('q', data=full_range_annular.mod_q)
    group.create_dataset('intensity', data=full_range_annular.intensity)

    for index, ring in enumerate(rings):
        group = debug_h5.create_group(f'slice {index}')
        group.create_dataset('q', data=ring.mod_q)
        group.create_dataset('intensity', data=ring.intensity)

    # close
    debug_h5.close()


def getWedgeSelection(data2d, q_min, q_delta, q_max, azimuthal_delta, peak_width=0.25, background_width=1.5,
                      signal_to_noise_min=2.):
    '''
    Calculate azimuthal binning ranges automatically based on finding peaks in the annular ring. The
    output of this is intended to be used in :py:func:`~drtsans.iq.select_i_of_q_by_wedge`.

    Parameters
    ==========
    data2d: ~drtsans.dataobjects.Azimuthal
    q_min: float
        The left bin boundary for the first Q-bin
    q_delta: float
        The size of the bins in Q
    q_max: float
        The left bin  boundary for the last Q-bin
    azimuthal_delta: float
        The size of the bins in azimuthal angle
    peak_width: float
        Percent of full-width-half-max (FWHM) of the peak to define the signal to be within when
        determining the final range for azimuthal binning.
    background_width: float
        Percent of full-width-half-max (FWHM) of the peak to define the background between peaks
        to be within when determining the final range for azimuthal binning.
    signal_to_noise_min: float
        Minimum signal to noise ratio for the data to be considered "fittable"

    Results
    =======
    ~list
      list containing 2 2-tuples as ``[((peak1_min, peak1_max), (peak2_min, peak2_max)), ((..., ...), (..., ...))]``
    '''
    # Bin azimuthal
    q, azimuthal_rings = _binInQAndAzimuthal(data2d, q_min=q_min, q_max=q_max, q_delta=q_delta,
                                             azimuthal_delta=azimuthal_delta)

    center_vec, fwhm_vec = _fitQAndAzimuthal(azimuthal_rings, q_bins=q,
                                             signal_to_noise_min=signal_to_noise_min, azimuthal_start=110.,
                                             maxchisq=1000.)

    # verify that the results didn't predict wedges larger than half of the data
    if np.any(np.array(fwhm_vec) > 360./2):
        values = ['{:.1f}deg'.format(value) for value in fwhm_vec]
        raise RuntimeError('Encountered large fwhm values: {}'.format(', '.join(values)))

    # convert to min and max ranges
    min_vec, max_vec = [], []

    min_vec.append(center_vec[0] - peak_width * fwhm_vec[0])
    max_vec.append(center_vec[0] + peak_width * fwhm_vec[0])

    min_vec.append(center_vec[0] + background_width * fwhm_vec[0])
    max_vec.append(center_vec[1] - background_width * fwhm_vec[1])

    min_vec.append(center_vec[1] - peak_width * fwhm_vec[1])
    max_vec.append(center_vec[1] + peak_width * fwhm_vec[1])

    min_vec.append(center_vec[1] + background_width * fwhm_vec[1])
    max_vec.append(center_vec[0] - background_width * fwhm_vec[0])

    # clean up the data to be in the form expected by select_i_of_q_by_wedge
    min_vec = np.array(min_vec)
    max_vec = np.array(max_vec)

    min_vec[min_vec < -90.] += 360.
    max_vec[max_vec < -90.] += 360.

    min_vec[min_vec > 270.] -= 360.
    max_vec[max_vec > 270.] -= 360.

    # put wedges on opposite sides together
    raw_wedges = list(zip(min_vec, max_vec))
    summing_wedges = []
    for i in range(len(raw_wedges) // 2):  # iterate over half the list
        summing_wedges.append((raw_wedges[i], raw_wedges[i+2]))

    return summing_wedges


def _binInQAndAzimuthal(data, q_min, q_delta, q_max, azimuthal_delta):
    '''This function bins the data in Qmod and azimuthal accoring to the supplied parameters. The maximum
    azimuthal is 540deg to allow for finding a peak at/near azimuthal=0deg.

    Parameters
    ==========
    data: ~drtsans.dataobjects.IQazimuthal
    q_min: float
        The left bin boundary for the first Q-bin
    q_delta: float
        The size of the bins in Q
    q_max: float
        The left bin  boundary for the last Q-bin
    azimuthal_delta: float
        The size of the bins in azimuthal

    Results
    =======
    tuple
        Histogram of ```(intensity, error, azimuthal_bins, q_bins)```
    '''
    # the bonus two steps is to get the end-point in the array
    q_bins = np.arange(q_min, q_max + q_delta, q_delta, dtype=float)

    # create azimuthal binning BinningParams takes number of steps
    azimuthal_offset = 0.5 * azimuthal_delta
    azimuthal_binning = BinningParams(0. - azimuthal_offset, 360. - azimuthal_offset,
                                      bins=int(360. / azimuthal_delta))
    # create the I(azimuthal) for each q-ring
    data_of_q_rings = []
    for qmin_ring, qmax_ring in zip(q_bins[:-1], q_bins[1:]):
        # bin into I(azimuthal)
        I_azimuthal = bin_annular_into_q1d(data, azimuthal_binning, qmin_ring, qmax_ring,
                                           BinningMethod.NOWEIGHT)

        # Create a copy of the arrays with the 360->540deg region repeated
        # ignore - delta_mod_q wavelength
        mod_q_new = determine_1d_linear_bins(x_min=0., x_max=540.+azimuthal_delta,
                                             bins=1 + int(540. / azimuthal_delta)).centers
        num_orig_bins = I_azimuthal.mod_q.size
        num_repeated_bins = mod_q_new.size - num_orig_bins

        intensity_new = np.zeros(mod_q_new.size)
        intensity_new[:num_orig_bins] = I_azimuthal.intensity
        intensity_new[-1 * num_repeated_bins:] = I_azimuthal.intensity[:num_repeated_bins]

        error_new = np.zeros(mod_q_new.size)
        error_new[:num_orig_bins] = I_azimuthal.error
        error_new[-1 * num_repeated_bins:] = I_azimuthal.error[:num_repeated_bins]

        I_azimuthal = IQmod(intensity=intensity_new, error=error_new, mod_q=mod_q_new)

        # append to the list of spectra
        data_of_q_rings.append(I_azimuthal)

    # return intensity, error, azimuthal_bins, q_bins TODO REMOVE
    return q_bins, data_of_q_rings


def _estimatePeakParameters(intensity, azimuthal, azimuthal_start, window_half_width):
    '''Estimate the peak parameters by determining a window around a bright point in the data then using
    moment calculations to estimate the parameters for a Gaussian that approximates the actual peak.
    This is done to aid the fitting which does better with better starting values.

    Parameters
    ==========
    intensity: numpy.ndarray
        Array of intensities. This must not have nans in it.
    azimuthal: numpy.ndarray
        Array of azimuthal angles centers
    azimuthal_start: float
        Starting guess of peak center
    window_half_width: float
        The window used is this amount on either side of what is determined to be the peak center

    Results
    =======
    tuple
        ``(intensity, center, sigma)`` where intensity is the full height including background
    '''
    # Look for the highest point in a section of the data. This is an iterative approach that starts with a window
    # centered at `azimuthal_start`, the repeats until the maximum within the window doesn't move more than 1deg.
    azimuthal_new = azimuthal_start  # where to search around
    azimuthal_last = azimuthal_start  # last known value
    while True:
        # determine new windows staying at least 90.deg inside the edges
        window_min = np.max((azimuthal_new - window_half_width, azimuthal.min() + 90.))
        window_max = np.min((azimuthal_new + window_half_width, azimuthal.max() - 90.))

        # create a search window around azimuthal_new
        left_index = azimuthal.searchsorted(window_min, side='right')
        right_index = azimuthal.searchsorted(window_max, side='right')
        # the highest value in the window
        max_value = intensity[left_index:right_index].max()
        # where that is in the window
        max_index = np.where(intensity[left_index:right_index] == max_value)[0].max() + left_index
        # update values
        azimuthal_last = azimuthal_new
        azimuthal_new = azimuthal[max_index]

        # stop searching if the value hasn't changed by less than one degree
        if abs(azimuthal_new - azimuthal_last) < 1.:
            break

    # now use the first two moments of the data within the window to give an improved center position (first moment)
    # and width (derived from second moment)

    # the position of the center of the peak is the first momement of the data. "mean" can be thought of as the
    # center of mass of the peak in azimuthal angle.
    mean = np.sum(intensity[left_index: right_index] * azimuthal[left_index: right_index]) \
        / np.sum(intensity[left_index: right_index])

    # the fit uses sigma rather than fwhm
    # calculate the second moment about the mean as an approximation to a Gaussian's "sigma" parameter
    sigma = np.sum(intensity[left_index: right_index] * np.square(azimuthal[left_index: right_index] - mean)) \
        / np.sum(intensity[left_index: right_index])
    sigma = np.sqrt(sigma)

    return max_value, mean, sigma


def _fitSpectrum(spectrum, q_value, signal_to_noise_min, azimuthal_start):
    '''Extract the peak fit parameters for the data. This is done by observing where 2 maxima are in the
    spectrum then fitting for the peak parameters. This makes the assumption that the two peaks are 180deg
    apart.

    Parameters
    ==========
    intensity: numpy.ndarray
        Array of intensities for the q_value
    error: numpy.ndarray
        Array of uncertainties for the q_value
    azimuthal_bins: numpy.ndarray
        Array of azimuthal angles bin boundaries
    q_value_start: float
        Label for Q value being used. This is used for improving error messages
    signal_to_noise_min: float
        Minimum signal to noise ratio for the data to be considered "fittable"
    azimuthal_start: float
        First position to look for peaks around

    Results
    =======
    dict
        dict[name] = (value, error) where all of the fit parameters are converted.
        f0 is background, then f1...fn are the fitted peaks
    '''
    # define a default window size based on the number of peaks the function supports
    # currently only two peaks that are approximately 180deg apart is supported
    NUM_PEAK = 2
    WINDOW_SIZE = 0.6 * (360. / NUM_PEAK)

    # filter out the nans
    mask = np.logical_not(np.isnan(spectrum.intensity))
    if np.sum(mask) < 10:  # do not allow fitting less points than there are parameters
        raise RuntimeError('Less than 8 points being fit with 7 parameters (found {} points)'.format(np.sum(mask)))

    # first estimate background as minimum value
    # this will be subtracted off from found intensities during estimation
    background = spectrum.intensity[mask].min()

    # check if there is signal to noise greater than 2
    # this calculation assumes that the background is positive
    signal_to_noise = np.sum(spectrum.intensity[mask]) / (float(np.sum(mask)) * background)
    if signal_to_noise < signal_to_noise_min:
        raise RuntimeError('Estimated signal to noise is smaller than {}: found {:.2f}'.format(signal_to_noise_min,
                                                                                               signal_to_noise))

    # start of what will eventually be the fit function by specifying the background
    function = ['name=FlatBackground,A0={}'.format(background)]

    # template for describing initial peak guess
    gaussian_str = 'name=Gaussian,Height={},PeakCentre={},Sigma={}'

    # guess where one peak might be, start with a window of WINDOW_SIZE each side around 110
    intensity_peak, azimuthal_first, sigma = _estimatePeakParameters(spectrum.intensity[mask],
                                                                     spectrum.mod_q[mask],
                                                                     azimuthal_start=azimuthal_start,
                                                                     window_half_width=WINDOW_SIZE)
    function.append(gaussian_str.format(intensity_peak-background, azimuthal_first, sigma))

    # assume the other peak is 360 / NUM_PEAK degrees away
    azimuthal_start = azimuthal_first + (360. / NUM_PEAK)
    intensity_peak, azimuthal_second, sigma = _estimatePeakParameters(spectrum.intensity[mask],
                                                                      spectrum.mod_q[mask],
                                                                      azimuthal_start=azimuthal_start,
                                                                      window_half_width=WINDOW_SIZE)
    function.append(gaussian_str.format(intensity_peak-background, azimuthal_second, sigma))

    # create workspace version of data
    # this includes the nans so `Fit` has to be told to ignore them
    q_azimuthal_workspace = spectrum.to_workspace()

    # fit the positions of the two suspected peaks
    fit_workspace_prefix = unique_workspace_dundername()
    try:
        fitresult = Fit(Function=';'.join(function), InputWorkspace=q_azimuthal_workspace, Output=fit_workspace_prefix,
                        StartX=spectrum.mod_q.min() + 90., EndX=spectrum.mod_q.min() + 90. + 360.,
                        OutputParametersOnly=True, IgnoreInvalidData=True)
    except RuntimeError as e:
        raise RuntimeError('Failed to fit Q={}'.format(q_value)) from e
    finally:
        DeleteWorkspace(q_azimuthal_workspace)

    if fitresult.OutputStatus != 'success':
        raise RuntimeError('Failed to fit Q={}'.format(q_value))

    # convert the table into a dict[name] = (value, error)
    result = {}
    for i in range(fitresult.OutputParameters.rowCount()):
        row = fitresult.OutputParameters.row(i)
        name = row['Name']
        if name.startswith('Cost function'):
            name = 'chisq'
        result[name] = (row['Value'], row['Error'])

    # delete fit results
    for label in ['_Parameters', '_NormalisedCovarianceMatrix']:
        DeleteWorkspace(Workspace=fit_workspace_prefix + label)

    return result


def _toPositionAndFWHM(fitresult, peak_label, maxchisq):
    '''Returns ((center, center_error), (width, width_error))
    If chisq is too large or any of the errors is nan, all return values are set to nan

    This also generates the weights as height / parameter_error'''
    if fitresult['chisq'][0] > maxchisq:
        return (np.nan, np.nan), (np.nan, np.nan)
    else:
        # height = fitresult[peak_label + '.Height']
        center = fitresult[peak_label + '.PeakCentre']
        fwhm = tuple([value * SIGMA_TO_FWHM for value in fitresult[peak_label + '.Sigma']])

    # Anything being nan suggests that a fit failed. Set everything to nan so they do not
    # contribute to the weighted average.
    if np.isnan(center[1]) or np.isnan(fwhm[1]) or center[1] == 0. or fwhm[1] == 0.:
        return (np.nan, np.nan), (np.nan, np.nan)

    # Weights for are height divided by uncertainty. This results in stronger peaks with lower fitting
    # uncertainty contributing more to the parameters in azimuthal angle.
    center = center[0], 1. / center[1]
    fwhm = fwhm[0], 1. / fwhm[1]

    return (center, fwhm)


def _weighted_position_and_width(peaks):
    '''For a set of peaks, calculate the weighted average position and weighted average fwhm

    Parameters
    ==========
    peaks: list
        [((position, weight), (fwhm, weight)), ...] Each is a peak for a single Q-value
    as a function of azimuthal angle.

    Results
    =======
    tuple
        (position, fwhm)
    '''
    if len(peaks) <= 0:
        raise RuntimeError('Encountered zero fitted peaks')
    pos_accum, pos_weight_accum = 0., 0.
    fwhm_accum, fwhm_weight_accum = 0., 0.
    for peak in peaks:
        # friendlier names
        pos, pos_weight = peak[0]  # position and weight
        fwhm, fwhm_weight = peak[1]  # fwhm and weight

        if np.isnan(pos_weight) or np.isnan(fwhm_weight) or pos_weight <= 0. or fwhm_weight <= 0.:
            continue  # don't use these points

        pos_accum += pos * pos_weight
        pos_weight_accum += pos_weight

        fwhm_accum += fwhm * fwhm_weight
        fwhm_weight_accum += fwhm_weight

    try:
        return (pos_accum / pos_weight_accum), (fwhm_accum / fwhm_weight_accum)
    except ZeroDivisionError as e:
        raise RuntimeError('Cannot determine fitted positions from zero weights') from e


def _fitQAndAzimuthal(azimuthal_rings, q_bins, signal_to_noise_min, azimuthal_start, maxchisq):
    '''Find the peaks in the azimuthal spectra, then combine them into
    two composite centers and fwhm. This is currently coded to only
    look for two peaks.

    Parameters
    ==========
    intensity: numpy.ndarray
        The intensity as a 2-d array of Q and azimuthal
    error: numpy.ndarray
        The uncertainties as a 2-d array of Q and azimuthal
    azimuthal_bins: numpy.ndarray
        Array of azimuthal angles bin boundaries
    q_bins: numpy.ndarray
        array of Q bin boundaries
    signal_to_noise_min: float
        Minimum signal to noise ratio for the data to be considered "fittable"
    azimuthal_start: float
        First position to look for peaks around
    maxchisq: float
        The maximum chisq value for a fit result to be used in calculating the composite peak

    Results
    =======
    list, list
        The first list is the peak centers, the second is the peak fwhm
    '''
    if len(azimuthal_rings) != len(q_bins) - 1:
        raise RuntimeError('must supply q-bin boundaries')

    # change to centers in Q for messages
    q_centers = 0.5 * (q_bins[:-1] + q_bins[1:])

    # select out a single spectrum
    peakResults = [[], []]
    q_centers_used = []

    index = 0
    used_index = list()
    unfit_message = ''

    for spectrum, q_center in zip(azimuthal_rings, q_centers):
        index += 1
        try:
            fitresult = _fitSpectrum(spectrum, q_center,
                                     signal_to_noise_min=signal_to_noise_min, azimuthal_start=azimuthal_start)
            newlyFittedPeaks = [_toPositionAndFWHM(fitresult, label, maxchisq) for label in ['f1', 'f2']]

            if np.isnan(newlyFittedPeaks[0][0][0]) or np.isnan(newlyFittedPeaks[1][0][0]):
                unfit_message += f'spectrum {index}: failed to fit peaks due to NaN in fit result\n'
                continue

            for i in range(len(peakResults)):
                peakResults[i].append(newlyFittedPeaks[i])
            q_centers_used.append(q_center)
            used_index.append(index - 1)
        except RuntimeError as e:
            unfit_message += 'spectrum {}: Not using information from Q-slice ({}A):'.format(index, q_center)
            unfit_message += f'Encountered runtime error: {e}\n'  # don't worry about it
            continue
        except ValueError as val_err:
            # in case user specifies a range containing no Q values
            unfit_message += f'Spectrum {index}: unable to fit peaks due to {val_err}\n'
            continue

    logger.notice(f'Q-rings used to determine overall wedges: {q_centers_used}')
    logger.information(f'used annular binning index: {used_index}')
    logger.notice(unfit_message)

    peakResults = [_weighted_position_and_width(peak) for peak in peakResults]

    # convert into parallel arrays of centers and fwhm
    center_list = []
    fwhm_list = []
    for center, fwhm in peakResults:
        center_list.append(center)
        fwhm_list.append(fwhm)

    logger.notice(f'Fitted peak centers: {center_list}\nFWHMs              : {fwhm_list}')

    return center_list, fwhm_list
