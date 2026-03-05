from typing import Tuple

from drtsans.polarization import polarized_sample


def resolve_slicing(reduction_input: dict) -> Tuple[bool, bool]:
    r"""
    Resolve if the reduction configuration parameters specify time or log slicing


    Parameters
    ----------
    reduction_input
        Dictionary of reduction configuration parameters

    Returns
    -------
    boolean values for time and log slicing, respectively

    Raises
    ------
    ValueError
        - If the sample input data is composed of more than one run
        - If both time and log slicing are ``True``
    """
    reduction_config = reduction_input["configuration"]
    timeslice, logslice = reduction_config["useTimeSlice"], reduction_config["useLogSlice"]
    if timeslice and logslice:
        raise ValueError("Can't do both time and log slicing")
    sample = reduction_input["sample"]["runNumber"]
    multiple_samples = len(sample.split(",")) > 1
    if (timeslice or logslice) and multiple_samples:
        raise ValueError("Can't do slicing on summed data sets")
    polarized = polarized_sample(reduction_input)
    if (timeslice and polarized) or (logslice and polarized):
        raise NotImplementedError("Time or log slicing on polarized data sets is not implemented yet")
    return timeslice, logslice, polarized
