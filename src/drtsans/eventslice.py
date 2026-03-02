from typing import Tuple


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
    parameters = reduction_input["configuration"]
    timeslice, logslice = parameters["useTimeSlice"], parameters["useLogSlice"]
    if timeslice and logslice:
        raise ValueError("Can't do both time and log slicing")
    if timeslice or logslice:
        sample = reduction_input["sample"]["runNumber"]
        if len(sample.split(",")) > 1:
            raise ValueError("Can't do slicing on summed data sets")
    return timeslice, logslice
