import enum
from mantid.kernel import ConfigService
from mantid.api import mtd

__all__ = ['InstrumentEnumName', 'instrument_enum_name', 'is_time_of_flight']


@enum.unique
class InstrumentEnumName(enum.Enum):
    r"""Unique names labelling each instrument"""
    UNDEFINED = None  # usually the dummy instrument used for testing
    BIOSANS = ConfigService.getFacility('HFIR').instrument('BIOSANS')
    EQSANS = ConfigService.getFacility('SNS').instrument('EQSANS')
    GPSANS = ConfigService.getFacility('HFIR').instrument('GPSANS')

    def __str__(self):
        return self.name


def instrument_enum_name(input_query):
    r"""
    Resolve the instrument name as a unique enumeration.

    Parameters
    ----------
    input_query: str,  ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        string representing a filepath, a valid instrument name, or a Mantid workspace containing an instrument

    Returns
    -------
    InstrumentEnumName
        The name of the instrument as one of the InstrumentName enumerations
    """
    string_to_enum = {'CG3': InstrumentEnumName.BIOSANS, 'BIOSANS': InstrumentEnumName.BIOSANS,
                      'EQ-SANS': InstrumentEnumName.EQSANS, 'EQSANS': InstrumentEnumName.EQSANS,
                      'CG2': InstrumentEnumName.GPSANS, 'GPSANS': InstrumentEnumName.GPSANS}
    # convert to a string
    name = str(input_query)

    if name in mtd:  # convert mantid workspace into a instrument string
        name = mtd[str(name)].getInstrument().getName()
    else:  # see if `name` contains any of the instrument labels
        for instrument_string_label in string_to_enum:
            if instrument_string_label in name.upper():
                name = instrument_string_label
                break

    return string_to_enum.get(name.upper(), InstrumentEnumName.UNDEFINED)


def is_time_of_flight(input_query):
    r"""
    Find if the instrument is a time-of-flight one

    Parameters
    ----------
    input_query: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventsWorkspace
        string representing a valid instrument name, or a Mantid workspace containing an instrument

    Returns
    -------
    bool
    """
    return instrument_enum_name(input_query) is InstrumentEnumName.EQSANS  # we only have one, for the moment
