from drtsans.mono.normalisation import (time, monitor)
# Functions exposed to the general user (public) API
__all__ = ['apply_solid_angle_correction', 'normalize']


def normalize(ws, normalization_type):
    """Normalize to time, monitor, or proton charge"""
    if normalization_type == "time":
        return time(input_workspace=ws)
    elif normalization_type == "monitor":
        return monitor(input_workspace=ws)
    else:
        raise NotImplementedError()