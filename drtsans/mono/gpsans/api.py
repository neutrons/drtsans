from drtsans import solid_angle_correction
from drtsans.mono.normalization import (normalize_by_time, normalize_by_monitor)

# Functions exposed to the general user (public) API
__all__ = ['apply_solid_angle_correction', 'normalize_by_flux']


def apply_solid_angle_correction(ws):
    """
        Apply solid angle correction
    """
    # correct implementation requires BOTH wing and tube detectors
    # see https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/107
    raise NotImplementedError()
    return solid_angle_correction(ws, detector_type='VerticalTube')
