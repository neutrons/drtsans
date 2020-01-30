import numpy as np
from drtsans.dataobjects import DataType, getDataType
from drtsans.settings import unique_workspace_dundername
# https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html
from mantid.simpleapi import CreateWorkspace
# https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
from mantid.simpleapi import DeleteWorkspace
# https://docs.mantidproject.org/nightly/algorithms/Fit-v1.html
from mantid.simpleapi import Fit

__all__ = ['getWedgeSelection']

# factor to convert from a sigma to full width half max
# https://docs.mantidproject.org/nightly/fitting/fitfunctions/Gaussian.html
SIGMA_TO_FWHM = 2. * np.sqrt(2. * np.log(2.))


def getWedgeSelection(data2d, q_min, q_delta, q_max, azimuthal_delta, peak_width=0.25, background_width=1.5):
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

    Results
    =======
    tuple
      tuple of tuples i.e. ``((angle_min1, angle_max1), (angle_min2, angle_max2), ...)``
    '''

    intensity, error, azimuthal, q = _binInQAndAzimuthal(data2d, q_min=q_min, q_max=q_max, q_delta=q_delta,
                                                         azimuthal_delta=azimuthal_delta)
    center_vec, fwhm_vec = _fitQAndAzimuthal(intensity, error, q_bins=q, azimuthal_bins=azimuthal)

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

    return list(zip(min_vec, max_vec))


def _toQmodAndAzimuthal(data):
    '''This function returns the values of qmod and azimuthal that are parallel
    to the original data array. It requiresthat the data is IQazimuthal

    Parameters
    ==========
    data: ~drtsans.dataobjects.Azimuthal

    Results
    =======
    tuple
        ```(qmod, azimuthal)``` with the same dimensionality as the data.intensity
        with Q in angstrom and azimuthal angle in degrees'''
    if not getDataType(data) == DataType.IQ_AZIMUTHAL:
        raise RuntimeError('Calculating qmod and azimuthal only works for IQazimuthal')

    # reshape the qx and qy if intensity array is 2d
    if len(data.intensity.shape) == 2 and len(data.qx.shape) == 1 and len(data.qy.shape) == 1:
        qx = np.tile(data.qx, (data.qy.shape[0], 1))
        qy = np.tile(data.qy, (data.qx.shape[0], 1)).transpose()
    else:
        qx = data.qx
        qy = data.qy

    # calculate q-scalar
    q = np.sqrt(np.square(qx) + np.square(qy))

    # azimuthal is expected to be positive so use cyclical nature of trig functions
    azimuthal = np.arctan2(qy, qx)
    azimuthal[azimuthal < 0.] += 2. * np.pi

    return q, np.rad2deg(azimuthal)


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
    q, azimuthal = _toQmodAndAzimuthal(data)

    # the bonus two steps is to get the end-point in the array
    q_bins = np.arange(q_min, q_max + q_delta, q_delta, dtype=float)
    # additional half-circle is to pick up things that are symmetric near azimuthal=0
    azimuthal_max = 540. + azimuthal_delta
    azimuthal_bins = np.arange(-.5 * azimuthal_delta, azimuthal_max, azimuthal_delta, dtype=float)

    # create output data
    intensity = np.zeros((azimuthal_bins.size-1, q_bins.size-1), dtype=float)
    error = np.zeros((azimuthal_bins.size-1, q_bins.size-1), dtype=float)

    # do the binning - twice around the circle
    # while this loops through the data twice, it does not require copying data
    for azimuthal_offset in [0., 360.]:  # first pass then second pass
        for q_val, azimuthal_val, i_val, e_val in zip(q.ravel(), azimuthal.ravel() + azimuthal_offset,
                                                      data.intensity.ravel(), data.error.ravel()):
            # stop searching past 540
            if azimuthal_val > azimuthal_max:
                break

            # find the correct bin in Q
            q_index = q_bins.searchsorted(q_val, side='right')
            if q_index >= q_bins.size or q_index == 0:
                print('FAILED TO FIND Q_BIN FOR {} from {}'.format(q_val, q_bins))
                continue

            # find the correct bin in azimuthal
            azimuthal_index = azimuthal_bins.searchsorted(azimuthal_val, side='right')
            if azimuthal_index >= azimuthal_bins.size or q_index == 0:
                print('FAILED TO FIND AZIMUTHAL_BIN FOR {}'.format(azimuthal_val))
                continue

            # increment the counts array
            intensity[azimuthal_index - 1, q_index - 1] += i_val
            error[azimuthal_index - 1, q_index - 1] += e_val

    # bins that didn't accumulate counts are nan
    mask = error == 0.
    intensity[mask] = np.nan
    error[mask] = np.nan

    return intensity, error, azimuthal_bins, q_bins


def _estimatePeakParameters(intensity, azimuthal, azimuthal_start, window_half_width):
    '''Estimate the peak parameters using statistical measures. This is done to aid the
    fitting which does better with better starting values.

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
    # find the maximum in the window and move the azimuthal_new to there
    # this assumes that the data maximum is close to the peak center
    azimuthal_new = azimuthal_start  # where to search around
    azimuthal_last = azimuthal_start  # last known value
    while True:
        # create a search window around azimuthal_new
        left_index = azimuthal.searchsorted(azimuthal_new - window_half_width, side='right')
        right_index = azimuthal.searchsorted(azimuthal_new + window_half_width, side='right')
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

    # now use simple statistics to estimate the peak parameters

    # the position of the center of the peak is the mean
    mean = np.sum(intensity[left_index: right_index] * azimuthal[left_index: right_index]) \
        / np.sum(intensity[left_index: right_index])

    # the fit uses sigma rather than fwhm
    sigma = np.sum(intensity[left_index: right_index] * np.square(azimuthal[left_index: right_index] - mean)) \
        / np.sum(intensity[left_index: right_index])
    sigma = np.sqrt(sigma)

    return max_value, mean, sigma


def _fitSpectrum(intensity, error, azimuthal_bins, q_value):
    '''Extract the peak fit parameters for the data. This is done by observing where 2 maxima are in the
    spectrum then fitting for the peak parameters

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

    Results
    =======
    dict
        dict[name] = (value, error) where all of the fit parameters are converted.
        f0 is background, then f1...fn are the fitted peaks
    '''
    # required signal to noise ratio
    SIGNAL_TO_NOISE_MIN = 2.0

    # centers are more useful for various things below
    azimuthal_centers = 0.5 * (azimuthal_bins[:-1] + azimuthal_bins[1:])

    # filter out the nans
    mask = np.logical_not(np.isnan(intensity))
    if np.sum(mask) < 10:
        raise RuntimeError('Less than 8 points being fit with 7 parameters (found {} points)'.format(np.sum(mask)))

    # first estimate background as minimum value
    # this will be subtracted off from found intensities during estimation
    background = intensity[mask].min()

    # check if there is signal to noise greater than 2
    # this calculation assumes that the background is positive
    signal_to_noise = np.sum(intensity[mask]) / (float(np.sum(mask)) * background)
    if signal_to_noise < SIGNAL_TO_NOISE_MIN:
        raise RuntimeError('Estimated signal to noise is smaller than {}: found {:.2f}'.format(SIGNAL_TO_NOISE_MIN,
                                                                                               signal_to_noise))

    # start of what will eventually be the fit function by specifying the background
    function = ['name=FlatBackground,A0={}'.format(background)]

    # template for describing initial peak guess
    gaussian_str = 'name=Gaussian,Height={},PeakCentre={},Sigma={}'

    # guess where one peak might be, start with a window of 110deg each side around 110
    intensity_peak, azimuthal_first, sigma = _estimatePeakParameters(intensity[mask], azimuthal_centers[mask],
                                                                     azimuthal_start=110., window_half_width=110.)
    function.append(gaussian_str.format(intensity_peak-background, azimuthal_first, sigma))

    # assume the other peak is 180 degrees away
    intensity_peak, azimuthal_second, sigma = _estimatePeakParameters(intensity[mask], azimuthal_centers[mask],
                                                                      azimuthal_start=azimuthal_first + 180.,
                                                                      window_half_width=110.)
    function.append(gaussian_str.format(intensity_peak-background, azimuthal_second, sigma))

    # create workspace version of data
    # this uses bin boundaries and includes the nans so `Fit` has to be told to ignore them
    q_azimuthal_workspace = unique_workspace_dundername()
    CreateWorkspace(DataX=azimuthal_bins, DataY=intensity, DataE=error, Nspec=1,
                    UnitX='Degrees', OutputWorkspace=q_azimuthal_workspace,
                    VerticalAxisUnit='MomentumTransfer', VerticalAxisValues=str(q_value),
                    Distribution=True, EnableLogging=False)

    # fit the positions of the two suspected peaks
    fit_workspace_prefix = unique_workspace_dundername()
    try:
        fitresult = Fit(Function=';'.join(function), InputWorkspace=q_azimuthal_workspace, Output=fit_workspace_prefix,
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
        height = fitresult[peak_label + '.Height']
        center = fitresult[peak_label + '.PeakCentre']
        fwhm = tuple([value * SIGMA_TO_FWHM for value in fitresult[peak_label + '.Sigma']])

    # it is bad if anything is nan
    if np.isnan(center[1]) or np.isnan(fwhm[1]) or center[1] == 0. or fwhm[1] == 0.:
        return (np.nan, np.nan), (np.nan, np.nan)

    # weights for are height divided by uncertainty
    center = center[0], height[0] / center[1]
    fwhm = fwhm[0], height[0] / fwhm[1]

    return (center, fwhm)


def _weighted_least_sq(peaks):
    '''For a set of peaks, calculate the weighted position and fwhm

    Parameters
    ==========
    peaks: list
        [((position, weight), (fwhm, weight)), ...]

    Results
    =======
    tuple
        (position, fwhm)
    '''
    pos_accum, pos_weight_accum = 0., 0.
    fwhm_accum, fwhm_weight_accum = 0., 0.
    for peak in peaks:
        # friendlier names
        pos, pos_w = peak[0]  # position and weight
        fwhm, fwhm_w = peak[1]  # fwhm and weight

        if np.isnan(pos_w) or np.isnan(fwhm_w) or pos_w == 0. or fwhm_w == 0.:
            continue  # don't use these points

        pos_accum += pos * pos_w
        pos_weight_accum += pos_w

        fwhm_accum += fwhm * fwhm_w
        fwhm_weight_accum += fwhm_w

    return (pos_accum / pos_weight_accum), (fwhm_accum / fwhm_weight_accum)


def _fitQAndAzimuthal(intensity, error, azimuthal_bins, q_bins, maxchisq=1000.):
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
    maxchisq: float
        The maximum chisq value for a fit result to be used in calculating the compositi peak

    Results
    =======
    list, list
        The first list is the peak centers, the second is the peak fwhm
    '''
    # change to centers in Q for messages
    q_centers = 0.5 * (q_bins[:-1] + q_bins[1:])

    # select out a single spectrum
    peakResults = [[], []]
    for spec_index, q_center in enumerate(q_centers):
        try:
            # print('Fitting spectrum {} with Q={}A'.format(spec_index, q_center))
            fitresult = _fitSpectrum(intensity.T[spec_index], error.T[spec_index], azimuthal_bins, q_center)
            newlyFittedPeaks = [_toPositionAndFWHM(fitresult, label, maxchisq) for label in ['f1', 'f2']]
            # TODO loop over the index
            if np.isnan(newlyFittedPeaks[0][0][0]) or np.isnan(newlyFittedPeaks[1][0][0]):
                continue

            for i in range(len(peakResults)):
                peakResults[i].append(newlyFittedPeaks[i])
        except RuntimeError as e:
            print('Not using information from Q-slice ({}A):'.format(q_center),
                  'Encountered runtime error:', e)  # don't worry about it
            continue

    peakResults = [_weighted_least_sq(peak) for peak in peakResults]

    # convert into parallel arrays of centers and fwhm
    center_list = []
    fwhm_list = []
    for center, fwhm in peakResults:
        center_list.append(center)
        fwhm_list.append(fwhm)

    return center_list, fwhm_list
