# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/mono/convert_to_q.py
import drtsans.mono.convert_to_q


def convert_to_q(ws, mode, **kwargs):
    r"""
    Convert a workspace with units of wavelength into a
    series of arrays: intensity, error, q (or q components),
    delta q (or delta q components), and wavelength.

    See https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/mono/convert_to_q.py
    for description
    """
    return drtsans.mono.convert_to_q.convert_to_q(ws, mode, **kwargs)
