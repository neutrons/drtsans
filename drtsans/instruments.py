import enum
from mantid.kernel import ConfigService
from mantid.api import mtd

__all__ = ['InstrumentEnumName', 'instrument_enum_name']


@enum.unique
class InstrumentEnumName(enum.Enum):
    r"""Unique names labelling each instrument"""
    BIOSANS = ConfigService.getFacility('HFIR').instrument('BIOSANS')
    EQSANS = ConfigService.getFacility('SNS').instrument('EQSANS')
    GPSANS = ConfigService.getFacility('HFIR').instrument('GPSANS')

    def __str__(self):
        return self.name


def instrument_enum_name(label):
    r"""
    Resolve the instrument name as a unique enumeration.

    Parameters
    ----------
    label: str, Workspace
        string representing a valid instrument name, or a Mantid workspace containing an instrument

    Returns
    -------
    InstrumentEnumName
        The name of the instrument as one of the InstrumentName enumerations
    """
    string_to_enum = {'CG3': InstrumentEnumName.BIOSANS, 'BIOSANS': InstrumentEnumName.BIOSANS,
                      'EQ-SANS': InstrumentEnumName.EQSANS, 'EQSANS': InstrumentEnumName.EQSANS,
                      'CG2': InstrumentEnumName.GPSANS, 'GPSANS': InstrumentEnumName.GPSANS}
    # convert to a string
    name = str(label)

    # convert mantid workspaces into a instrument string
    if name in mtd:
        name = mtd[str(name)].getInstrument().getName()

    # dict only checks for uppercase names
    name = name.upper()

    # We want the enum representation of an instrument name
    if name in string_to_enum.keys():
        return string_to_enum[name]
    else:
        raise ValueError('Do not know how to convert "{}" to InstrumentName'.format(label))
