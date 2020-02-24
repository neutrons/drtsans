import numpy as np

# Functions exposed to the general user (public) API
__all__ = ['attenuation_factor']


def attenuation_factor(workspace):
    """This get the wavelength and attenuator value from the workspace
    logs then calculates the attenuation factor based on the fitted
    parameters for each different attenuator based on the equation
    scale = A * exp(-B * λ) + C

    The attenuation scale factor and the error for this is returned.

    The attenuator value pulled from the logs is mapped to the attenuator name by:

    0: "Undefined"
    1: "Close"
    2: "Open",
    3: "x3",
    4: "x30",
    5: "x300",
    6: "x2k",
    7: "x10k",
    8: "x100k"

    If the attenuator is one of Undefined, Close or Open then a scale factor of 1 with error 0 is returned

    """

    # Get attenuator value
    attenuator = workspace.getRun()['attenuator'].firstValue()

    if attenuator < 3:  # Undefined, Close, Open
        # return scale factor of 1 and error 0
        return 1, 0

    # The fitted attenuator parameters for the equation A * exp(-B * λ) + C
    # Provided by Lisa Debeer-Schmitt, 2020-02-21
    # In the following format (Amp, Amp Err, exp const, exp const err,  bkgd, bkgd err)
    attenuators = {3:  # x3
                   (0.3733459538730628, 0.008609717163831113,
                    0.08056544906925872, 0.008241433507695071,
                    0.0724341919138054, 0.01125160779959418),
                   4:  # x30
                   (0.11696573514650677, 0.006304060228295941,
                    0.25014801934427583, 0.01012469612884642,
                    0.003696051816711061, 0.0003197928933191539),
                   5:  # x300
                   (0.03280144048179223, 0.016878797252331778,
                    0.41335832054829935, 0.08308901706841165,
                    0.00020409078915193703, 6.186822388235207e-05),
                   6:  # x2k
                   (0.016574500241863355, 0.005593117048321498,
                    0.598755990591887, 0.06477682452109962,
                    7.602576960518972e-05, 1.3575163517990533e-05),
                   7:  # x10k
                   (0.00563013075327734, 0.0005203265715819975,
                    0.6961581698084675, 0.01938010154115584,
                    1.3123049075167468e-05, 1.5266828654446554e-06),
                   8:  # x100k
                   (0.1439135754790426, 0.005573924841205431,
                    0.30824770207752383, 0.011967728090637404,
                    0.006739099792400909, 0.0007076026868930973)}

    # Get wavelength from workspace
    λ = workspace.getRun()['wavelength'].timeAverageValue()

    # Call the attenuation function and return results
    return _attenuation_factor(*attenuators[attenuator], λ)


def _attenuation_factor(A, A_e, B, B_e, C, C_e, λ):
    """
    This calculates the function
        A * exp(-B * λ) + C
    along with the uncertainty
    """
    scale = A * np.exp(-B * λ) + C
    scale_error_Amp = np.exp(-B * λ)
    scale_error_exp_const = A * np.exp(-B * λ) * (-λ)
    scale_error_bkgd = 1
    scale_error = np.sqrt((scale_error_Amp * A_e)**2 + (scale_error_exp_const * B_e)**2 + (scale_error_bkgd * C_e)**2)
    return scale, scale_error
