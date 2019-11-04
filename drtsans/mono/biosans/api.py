from drtsans.mono.normalization import (normalize_by_time, normalize_by_monitor)

# Functions exposed to the general user (public) API
__all__ = ['normalize_by_flux']


def normalize_by_flux(ws, method):
    """Normalize to time or monitor"""
    if method == "time":
        return normalize_by_time(input_workspace=ws)
    elif method == "monitor":
        return normalize_by_monitor(input_workspace=ws)
    else:
        raise NotImplementedError()
