import numpy as np


def mask_as_numpy_array(w, invert=False):
    """Return mask in pixels as numpy array of bool. Items are True if masked
    :param w: input workspace
    :param invert: invert the array, items are True if unmasked
    :return: numpy.ndarray(bool)
    """
    mask = [w.getDetector(i).isMasked()
            for i in range(w.getNumberHistograms())]
    mask = np.asarray(mask)
    return mask if invert is False else np.invert(mask)


def masked_indexes(w, invert=False):
    """List of masked workspaces indexes
    :param w: input workspace
    :param invert: Return list of unmasked workspace indexes if True
    :return: numpy.ndarray(bool)
    """
    mask = mask_as_numpy_array(w, invert=invert)
    return np.where(mask)[0]
