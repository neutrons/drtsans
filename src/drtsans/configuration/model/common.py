from enum import Enum
from typing import Annotated, List, Optional, Union

from pydantic import BaseModel, Field, PositiveFloat, PositiveInt
from pydantic.dataclasses import dataclass
from pydantic.functional_validators import AfterValidator

from drtsans.configuration.model.validators import (
    validate_run_number_type,
    validate_safe_string_float,
    validate_safe_string_positive_float,
    validate_transmission_type,
)

### Definitions and Type Aliases

RunNumberTypes = Annotated[
    Union[str, PositiveInt, List[Union[str, PositiveInt]]], AfterValidator(validate_run_number_type)
]
RunNumberOptionalTypes = Optional[Union[RunNumberTypes, List]]
TransmissionValueTypes = Annotated[Union[str, float, None], AfterValidator(validate_transmission_type)]
SafeString = Optional[str]
SafeStringFloat = Annotated[Union[str, float, None], AfterValidator(validate_safe_string_float)]
SafeStringPositiveFloat = Annotated[
    Union[str, PositiveFloat, None], AfterValidator(validate_safe_string_positive_float)
]


# class EventsLoaderOptions:
#     """Options that can be passed to Mantid algorithms LoadEventNexus or LoadEvenAsWorkspace2D

#     Examples:
#     ---------

#     {'FilterByTofMin': 1000.0, 'FilterByTofMax': 16000.0}
#     """

#     def __init__(self, options: dict = {}) -> None:
#         self.__dict__ = options

#     def __str__(self):
#         return f"EventsLoaderOptions: {self.__dict__}"


@dataclass
class Transmission:
    """The transmission for the sample, run number or value (0 < value <=1). Can be empty

    Attributes:
    -----------
    run_number: RunNumberOptionalTypes
        The run number(s) for the transmission. This can be a single run number, a list of run numbers,
        or a comma-separated string of run numbers.
    value: TransmissionValueTypes
        The transmission value (0 < value <=1)
    error_tolerance: TransmissionValueTypes
        Maximum allowed relative error in the calculated transmission (0 < value <=1)
    """

    run_number: RunNumberOptionalTypes = Field(RunNumberOptionalTypes, alias="runNumber")
    value: TransmissionValueTypes = Field(TransmissionValueTypes, alias="value")
    error_tolerance: TransmissionValueTypes = Field(None, alias="errorTolerance")


@dataclass
class LMFitParameter:
    """Fit parameters

    Attributes:
    -----------
    value: float
        Numeric value of the parameter
    vary: bool
        Whether the parameter should be allowed to vary during the fit
    min: float
        Lower bound for the parameter
    max: float
        Upper bound for the parameter
    expr: str
        Mathematical expression used to constrain the parameter during the fit
    """

    value: float
    vary: bool
    min: float
    max: float
    expr: str


@dataclass
class GaussianCenteringOptions:
    """Gaussian centering options

    Attributes:
    -----------
    amp: Optional[LMFitParameter]
        Amplitude of the Gaussian function. Default: ws.extractY().max()
    sigma_x: Optional[LMFitParameter]
        X spead of the Gaussian function. Default: 0.01
    sigma_y: Optional[LMFitParameter]
        Y spead of the Gaussian function. Default: 0.01
    theta: Optional[LMFitParameter]
        Clockwise rotation angle of Gaussian function. Default: 0
    center_x: Optional[LMFitParameter]
        Estimate for the beam center in X [m]. Default: 0
    center_y: Optional[LMFitParameter]
        Estimate for the beam center in Y [m]. Default: 0
    """

    amp: Optional[LMFitParameter] = None
    sigma_x: Optional[LMFitParameter] = 0.01
    sigma_y: Optional[LMFitParameter] = 0.01
    theta: Optional[LMFitParameter] = 0.0
    center_x: Optional[LMFitParameter] = 0.0
    center_y: Optional[LMFitParameter] = 0.0


@dataclass
class COMCenteringOptions:
    """Center of mass centering options

    Attributes:
    -----------
    center_x: float
        Estimate for the beam center in X [m] (default: 0)
    center_y: float
        Estimate for the beam center in Y [m] (default: 0)
    tolerance: SafeStringPositiveFloat
        Tolerance for center of mass position between each iteration (default: 0.00125)
    direct_beam: bool
        If true, a direct beam calculation will be performed,
        otherwise, the center of mass of the scattering data will be computed by excluding the beam area
    beam_radius: SafeStringPositiveFloat
        Radius of the beam area [m], used to exclude the beam when
        calculating the center of mass of the scattering pattern
    integration_radius: Optional[SafeStringPositiveFloat]
        Radius in meters, used for including the pixels that are within the area defined by the integration_radius,
        when calculating the center of mass for asymmetric detector arrays
    """

    center_x: float = 0.0
    center_y: float = 0.0
    direct_beam: Optional[bool] = None
    beam_radius: Optional[SafeStringPositiveFloat] = None
    tolerance: SafeStringPositiveFloat = 0.00125
    integration_radius: Optional[SafeStringPositiveFloat] = None


### Models


class InstrumentName(str, Enum):
    """Instrument names"""

    BIOSANS = "BIOSANS"
    EQSANS = "EQSANS"
    GPSANS = "GPSANS"


inst = InstrumentName("BIOSANS")


class CenteringMethod(str, Enum):
    """Centering methods"""

    GAUSSIAN = "GAUSSIAN"
    CENTER_OF_MASS = "CENTER_OF_MASS"


class Sample(BaseModel):
    """Sample information

    Attributes:
    -----------
    run_number: RunNumberTypes
        The run number(s) for the sample. This can be a single run number, a list of run numbers,
        or a comma-separated string of run numbers.
    load_options: Optional[EventsLoaderOptions]
        Options that can be passed to Mantid algorithms LoadEventNexus or LoadEvenAsWorkspace2D.
        Consult the documentation of these algorithms to find the available options.
        These options will take effect when loading the sample run(s), but not other runs
    thickness: Union[str, PositiveFloat]
        Sample thickness in cm
    transmission: Transmission
        The transmission for the sample
    """

    runNumber: RunNumberTypes
    loadOptions: Optional[dict] = {}
    thickness: Union[str, PositiveFloat]  # minLength = 1
    transmission: Transmission


class Background(BaseModel):
    """Background information

    Attributes:
    -----------
    run_number: RunNumberTypes
        The run number(s) for the background. This can be a single run number, a list of run numbers,
        or a comma-separated string of run numbers.
    transmission: Transmission
        The transmission for the background
    """

    runNumber: RunNumberOptionalTypes
    transmission: Transmission


class BeamCenter(BaseModel):
    """Beam center information

    Attributes:
    -----------
    run_number: RunNumberOptionalTypes

    """

    runNumber: RunNumberOptionalTypes
    method: Optional[CenteringMethod] = None
    gaussian_centering_options: Optional[GaussianCenteringOptions] = None
    com_centering_options: Optional[COMCenteringOptions] = None


class ReductionParameters(BaseModel):
    """Reduction parameters for DRTSANS data reduction

    Attributes:
    -----------
    instrument: InstrumentName
        The instrument used for the data collection
    sample: Sample
        Sample information
    background: Background
        Background information
    beam_center: BeamCenter
        Beam center information
    output_name: str
        The prefix for all output filenames, such as {outputFileName}_Iqxqy.png. It cannot be left empty.
    """

    # instrument: InstrumentName
    instrument: str
    sample: Sample
    background: Background
    beam_center: BeamCenter
    output_name: str
