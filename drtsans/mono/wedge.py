import numpy as np
from drtsans.dataobjects import DataType, getDataType
from drtsans.settings import unique_workspace_dundername
# https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html
from mantid.simpleapi import CreateWorkspace
# https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
from mantid.simpleapi import DeleteWorkspace
# https://docs.mantidproject.org/nightly/algorithms/Fit-v1.html
from mantid.simpleapi import Fit

# factor to convert from a sigma to full width half max
# https://docs.mantidproject.org/nightly/fitting/fitfunctions/Gaussian.html
SIGMA_TO_FWHM = 2. * np.sqrt(2. * np.log(2.))


def _toQmodAndPhi(data):
    '''This function returns the values of qmod and phi that are parallel
    to the original data array. It assumes that the data is IQazimuthal

    Parameters
    ==========
    data: ~drtsans.dataobjects.Azimuthal

    Results
    =======
    tuple
        ```(qmod, phi)``` with the same dimensionality as the data.intensity
        with Q in angstrom and phi in degrees'''
    if not getDataType(data) == DataType.IQ_AZIMUTHAL:
        raise RuntimeError('Calculating qmod and phi only works for IQazimuthal')

    # reshape the qx and qy if intensity array is 2d
    if len(data.intensity.shape) == 2 and len(data.qx.shape) == 1 and len(data.qy.shape) == 1:
        qx = np.tile(data.qx, (data.qy.shape[0], 1))
        qy = np.tile(data.qy, (data.qx.shape[0], 1)).transpose()
    else:
        qx = data.qx
        qy = data.qy

    # calculate q-scalar
    q = np.sqrt(np.square(qx) + np.square(qy))

    # phi is expected to be positive so use cyclical nature of trig functions
    phi = np.arctan2(qy, qx)
    phi[phi < 0.] += 2. * np.pi

    return q, np.rad2deg(phi)


def _binInQAndPhi(data, q_min, q_delta, q_max, phi_delta):
    '''This function bins the data in Qmod and phi accoring to the supplied parameters. The maximum phi is
    540deg to allow for finding a peak at/near phi=0deg.f

    Parameters
    ==========
    data: ~drtsans.dataobjects.Azimuthal
    q_min: float
        The left bin boundary for the first Q-bin
    q_delta: float
        The size of the bins in Q
    q_max: float
        The left bin  boundary for the last Q-bin
    phi_delta: float
        The size of the bins in phi

    Results
    =======
    tuple
        Histogram of ```(intensity, error, phi_bins, q_bins)```
    '''
    q, phi = _toQmodAndPhi(data)

    # the bonus two steps is to get the end-point in the array
    q_bins = np.arange(q_min, q_max + q_delta, q_delta, dtype=float)
    # additional half-circle is to pick up things that are symmetric near phi=0
    phi_max = 540. + phi_delta
    phi_bins = np.arange(-.5 * phi_delta, phi_max, phi_delta, dtype=float)

    # create output data
    intensity = np.zeros((phi_bins.size-1, q_bins.size-1), dtype=float)
    error = np.zeros((phi_bins.size-1, q_bins.size-1), dtype=float)

    # do the binning - twice around the circle
    # while this loops through the data twice, it does not require copying data
    for phi_offset in [0., 360.]:  # first pass then second pass
        for q_val, phi_val, i_val, e_val in zip(q.ravel(), phi.ravel() + phi_offset,
                                                data.intensity.ravel(), data.error.ravel()):
            # stop searching past 540
            if phi_val > phi_max:
                break

            # find the correct bin in Q
            q_index = q_bins.searchsorted(q_val, side='right')
            if q_index >= q_bins.size or q_index == 0:
                print('FAILED TO FIND Q_BIN FOR {} from {}'.format(q_val, q_bins))
                continue

            # find the correct bin in phi
            phi_index = phi_bins.searchsorted(phi_val, side='right')
            if phi_index >= phi_bins.size or q_index == 0:
                print('FAILED TO FIND PHI_BIN FOR {}'.format(phi_val))
                continue

            # increment the counts array
            intensity[phi_index - 1, q_index - 1] += i_val
            error[phi_index - 1, q_index - 1] += e_val

    # bins that didn't accumulate counts are nan
    mask = error == 0.
    intensity[mask] = np.nan
    error[mask] = np.nan

    return intensity, error, phi_bins, q_bins


def _estimatePeakParameters(intensity, phi, phi_start, window_half_width):
    '''Estimate the peak parameters using statistical measures. This is done to aid the
    fitting which will do better with better starting values.

    Parameters
    ==========
    intensity: numpy.ndarray
        Array of intensities. This must not have nans in it.
    phi: numpy.ndarray
        Array of phi angles centers
    phi_start: float
        Starting guess of peak center
    window_half_width: float
        The window used is this amount on either side of what is determined to be the peak center

    Results
    =======
    tuple
        ``(intensity, center, sigma)`` where intensity is the full height including background
    '''
    # find the maximum in the window and move the phi_new to there
    # this assumes that the data maximum is close to the peak center
    phi_new = phi_start  # where to search around
    phi_last = phi_start  # last known value
    while True:
        # create a search window around phi_new
        left_index = phi.searchsorted(phi_new - window_half_width, side='right')
        right_index = phi.searchsorted(phi_new + window_half_width, side='right')
        # the highest value in the window
        max_value = intensity[left_index:right_index].max()
        # where that is in the window
        max_index = np.where(intensity[left_index:right_index] == max_value)[0].max() + left_index
        # update values
        phi_last = phi_new
        phi_new = phi[max_index]

        # stop searching if the value hasn't changed by less than one degree
        if abs(phi_new - phi_last) < 1.:
            break

    # now use simple statistics to estimate the peak parameters

    # the position of the center of the peak is the mean
    mean = np.sum(intensity[left_index: right_index] * phi[left_index: right_index]) \
        / np.sum(intensity[left_index: right_index])

    # the fit uses sigma rather than fwhm
    sigma = np.sum(intensity[left_index: right_index] * np.square(phi[left_index: right_index] - mean)) \
        / np.sum(intensity[left_index: right_index])
    sigma = np.sqrt(sigma)

    return max_value, mean, sigma


def _fitSpectrum(intensity, error, phi_bins, q_value):
    '''Extract the peak fit parameters for the data. This is done by observing where 2 maxima are in the
    spectrum then fitting for the peak parameters

    Parameters
    ==========
    intensity: numpy.ndarray
        Array of intensities for the q_value
    error: numpy.ndarray
        Array of uncertainties for the q_value
    phi_bins: numpy.ndarray
        Array of phi angles bin boundaries
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
    phi_centers = 0.5 * (phi_bins[:-1] + phi_bins[1:])

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
    intensity_peak, phi_first, sigma = _estimatePeakParameters(intensity[mask], phi_centers[mask],
                                                               phi_start=110., window_half_width=110.)
    function.append(gaussian_str.format(intensity_peak-background, phi_first, sigma))

    # assume the other peak is 180 degrees away
    intensity_peak, phi_second, sigma = _estimatePeakParameters(intensity[mask], phi_centers[mask],
                                                                phi_start=phi_first + 180., window_half_width=110.)
    function.append(gaussian_str.format(intensity_peak-background, phi_second, sigma))

    # create workspace version of data
    # this uses bin boundaries and includes the nans so `Fit` has to be told to ignore them
    q_phi_workspace = unique_workspace_dundername()
    CreateWorkspace(DataX=phi_bins, DataY=intensity, DataE=error, Nspec=1,
                    UnitX='Degrees', OutputWorkspace=q_phi_workspace,
                    VerticalAxisUnit='MomentumTransfer', VerticalAxisValues=str(q_value),
                    Distribution=True, EnableLogging=False)

    # fit the positions of the two suspected peaks
    fit_workspace_prefix = unique_workspace_dundername()
    try:
        fitresult = Fit(Function=';'.join(function), InputWorkspace=q_phi_workspace, Output=fit_workspace_prefix,
                        OutputParametersOnly=True, IgnoreInvalidData=True)
    except RuntimeError as e:
        raise RuntimeError('Failed to fit Q={}'.format(q_value)) from e
    finally:
        DeleteWorkspace(q_phi_workspace)

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

    # weights for
    center = center[0], height[0] / center[1]
    fwhm = fwhm[0], height[0] / fwhm[1]

    return (center, fwhm)


def _weighted_least_sq(peaks):
    '''For a set of peaks, calculate the weighted positiojn and fwhm

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
        pos, pos_w = peak[0]  # position and error
        fwhm, fwhm_w = peak[1]  # fwhm and error

        if np.isnan(pos_w) or np.isnan(fwhm_w) or pos_w == 0. or fwhm_w == 0.:
            continue  # don't use these points

        pos_accum += pos * pos_w
        pos_weight_accum += pos_w

        fwhm_accum += fwhm * fwhm_w
        fwhm_weight_accum += fwhm_w

    return (pos_accum / pos_weight_accum), (fwhm_accum / fwhm_weight_accum)


def _fitQAndPhi(intensity, error, phi_bins, q_bins, maxchisq=1000.):
    # change to centers in Q for messages
    q_centers = 0.5 * (q_bins[:-1] + q_bins[1:])

    # select out a single spectrum
    peakResults = [[], []]
    for spec_index, q_center in enumerate(q_centers):
        try:
            # print('Fitting spectrum {} with Q={}A'.format(spec_index, q_center))
            fitresult = _fitSpectrum(intensity.T[spec_index], error.T[spec_index], phi_bins, q_center)
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

    return tuple(peakResults)
