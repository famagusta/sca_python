"""
Some helper functions that are used in the tutorial files
"""
import numpy as np


def truncate(alignments, frac_alpha_cutoff):
    """
    truncate alignments to sequence positions with percentage gap frequency no greater than
    100 * (1 - frac_alpha_cutoff)
    :param alignments:
    :param frac_alpha_cutoff:
    :return: alignments array after truncation
    """
    n_alignments = len(alignments)
    is_letter = np.core.defchararray.isalpha(alignments)
    frac_letters = is_letter.sum(0) * 1.0 / n_alignments
    return alignments[:, frac_letters >= frac_alpha_cutoff]
