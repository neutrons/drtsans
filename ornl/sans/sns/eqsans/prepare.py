from __future__ import (absolute_import, division, print_function)


from ornl.settings import load_run
from ornl.sans.sns.eqsans.beam_finder import direct_beam_center
from mantid.simpleapi import (SANSMaskDTP, FindCenterOfMassPosition)


def prepare_direct_beam_center(direct_beam, mask=None,
                               finder=FindCenterOfMassPosition,
                               finder_kwargs=None):
    r"""Recipe to find the beam center coordinates from a  direct beam run

    Current limitation: `mask` has to be a list of tube indexes

    Parameters
    ----------
    direct_beam: str
        Run number for the direct beam run
    mask: str,
        list of tubes (separated by comma) to mask.
    finder: function
        Method to find the beam center
    finder_kwargs: dict
        Additional options for the finder method
    """
    if finder_kwargs is None:
        finder_kwargs = {}
    w = load_run(direct_beam, 'EQSANS')
    if mask is not None:
        SANSMaskDTP(InputWorkspace=w, Tube=mask)
    return direct_beam_center(w, finder, finder_kwargs)


def prepare_dark_current():
    r"""Recipe to normalize the dark current"""
    pass


def prepare_sensitivity():
    r"""Recipe to generate a sensitivity file from a flood run"""
    pass


def prepare_scattering():
    r"""Prepare a run containing an element capable of scattering neutrons"""
    pass


def prepare_sample_scattering():
    r"""Prepare a run containing a sample"""
    return prepare_scattering()


def prepare_background_scattering():
    r"""Prepare a run containing an empty sample holder, or a background
    material"""
    return prepare_scattering()


def prepare_transmission():
    r"""
    Prepare a run containing an element capable of scattering neutrons
    attenuator plus an attenuator
    """
    pass


def prepare_sample_transmission():
    r"""Prepare a run containing a sample plus an attenuator"""
    return prepare_transmission()


def prepare_background_transmission():
    r"""
    Prepare a run containing an empty sample holder or a background
    material, plus an attenuator
    """
    return prepare_transmission()
