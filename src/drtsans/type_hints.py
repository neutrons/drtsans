"""
Compendium of custom type hints.

"""

# standard library imports
from collections.abc import Callable
from typing import TypeAlias, Union

# third party imports
import mantid


"""Any type of Mantid workspace, including its name"""
MantidWorkspace: TypeAlias = Union[str, mantid.api.Workspace]

"""Delayed moderator emission time as a function of wavelength, or no correction"""
EmissionDelay: TypeAlias = Callable[[float], float] | None
