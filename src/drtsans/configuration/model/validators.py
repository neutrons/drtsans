import re
from typing import Union, TypeVar

from pydantic import PositiveFloat


T = TypeVar("T")


def validate_safe_string_float(value: Union[str, float, None]) -> float:
    if value is None:
        return None
    if isinstance(value, str):
        if not re.match(r"^$|^[0-9]*.[0-9]*$", value):
            raise ValueError("Value must be a float")
    elif isinstance(value, float):
        return value
    else:
        raise ValueError("Value must be a float")
    return float(value)


def validate_safe_string_positive_float(value: Union[str, PositiveFloat, None]) -> float:
    if value is None:
        return None
    if isinstance(value, str):
        if not re.match(r"^$|^[0-9]*.[0-9]*$", value):
            raise ValueError("Value must be a float")
    elif isinstance(value, float):
        if value <= 0:
            raise ValueError("Value must be a positive float")
    else:
        raise ValueError("Value must be a float")
    return float(value)


def validate_transmission_type(value: Union[str, float, None]) -> float:
    ERROR = "Value must be a non-zero float between 0 and 1, or a string representation of a float"
    if value is None or value == "":
        return None
    if isinstance(value, str):
        try:
            value = float(value)
        except ValueError:
            raise ValueError(ERROR)
    elif isinstance(value, float):
        pass
    else:
        raise TypeError(ERROR)
    if value <= 0 or value > 1:
        raise ValueError(ERROR)
    return float(value)


def validate_run_number_type(value):
    ERROR = "Value must be a non-empty string, positive int, or non-empty list of strings or ints"
    if isinstance(value, str):
        if len(value) == 0:
            raise ValueError(ERROR)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, list):
        if len(value) == 0:
            raise ValueError(ERROR)
    return value
