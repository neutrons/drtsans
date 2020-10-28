import numpy as np
import lmfit


# defining 2D Gaussian fitting functions
def _Gaussian2D(x1, y1, amp, sigma_x, sigma_y, theta, x0, y0):
    a = np.cos(theta)**2/(2.*sigma_x**2) + np.sin(theta)**2/(2.*sigma_y**2)
    b = -np.sin(2.*theta)/(4.*sigma_x**2) + np.sin(2.*theta)/(4.*sigma_y**2)
    c = np.sin(theta)**2/(2.*sigma_x**2) + np.cos(theta)**2/(2.*sigma_y**2)
    amplitude = amp/(np.sqrt(2.*np.pi)*np.sqrt(sigma_x*sigma_y))
    val = amplitude * np.exp(-(a*(x1-x0)**2 + 2.*b*(x1-x0)*(y1-y0) + c*(y1-y0)**2))
    return val


def find_beam_center_gaussian(ws, parameters={}):
    # fitting 2D gaussian to center data
    ws_size = ws.getNumberHistograms()
    x = np.empty(ws_size)
    y = np.empty(ws_size)
    intes = np.empty(ws_size)
    intes_err = np.empty(ws_size)
    keep = np.empty(ws_size, dtype=np.bool_)

    for i, si in enumerate(ws.spectrumInfo()):
        pos = si.position
        x[i] = pos.X()
        y[i] = pos.Y()
        keep[i] = not si.isMasked and np.isfinite(ws.readY(i)[0])
        intes[i] = ws.readY(i)[0]
        intes_err[i] = ws.readE(i)[0]

    x = x[keep]
    y = y[keep]
    intes = intes[keep]
    intes_err = intes_err[keep]

    model = lmfit.Model(_Gaussian2D, independent_vars=["x1", "y1"],
                        param_names=["amp", "sigma_x", "sigma_y", "theta", "x0", "y0"])

    params = lmfit.Parameters()
    for key, value in parameters.items():
        params.add(key, **value)

    if 'amp' not in params:
        params.add('amp', value=ws.extractY().max())
    if 'sigma_x' not in params:
        params.add('sigma_x', value=0.01, min=np.finfo(float).eps)  # width in x
    if 'sigma_y' not in params:
        params.add('sigma_y', value=0.01, min=np.finfo(float).eps)  # width in y
    if 'theta' not in params:
        params.add('theta', value=0., min=-np.pi/2., max=np.pi/2.)
    if 'x0' not in params:
        params.add('x0', value=0.)
    if 'y0' not in params:
        params.add('y0', value=0.)
    results = model.fit(intes, x1=x, y1=y, weights=1./intes_err, params=params)
    return results.params['x0'].value, results.params['y0'].value
