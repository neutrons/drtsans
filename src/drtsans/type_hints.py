"""
Compendium of custom type hints.

"""

# standard library imports
from typing import TypeAlias, Union

# third party imports
import mantid


"""Any type of Mantid workspace, including its name"""
MantidWorkspace: TypeAlias = Union[str, mantid.api.Workspace]