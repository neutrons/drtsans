import pytest
import numpy as np


@pytest.mark.parametrize('workspace_with_instrument', [{'Nx': 3, 'Ny': 3}], indirect=True)
def test_apply_simple_sensitivity(workspace_with_instrument):
    """
    Testing section 5 in the master document
    Apply sensitivity to a 3x3 workspace. Check if the output is masked where sensitivity is masked
    Functions to test: drtsans.sensitivity.apply_sensitivity_correction
    Underlying Mantid algorithms:
        Divide https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html
        MaskDetectorsIf https://docs.mantidproject.org/nightly/algorithms/MaskDetectorsIf-v1.html
    See also https://docs.mantidproject.org/nightly/concepts/ErrorPropagation.html
    dev - Andrei Savici <saviciat@ornl.gov>
    SME - Venky Pingali <pingalis@ornl.gov>
    """
    data = np.array([[7., 8., 12.],
                     [10., 17., 13.],
                     [10., 10., 9.]])

    data_error = np.sqrt(data)

    # create workspaces
    data_ws = workspace_with_instrument(axis_values=[6.],  # fake wavelength
                                        intensities=data,
                                        uncertainties=data_error,
                                        view='pixel')

    assert data_ws is not None
