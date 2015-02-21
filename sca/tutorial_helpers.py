import numpy as np


def truncate(alignments, frac_alpha_cutoff):
    n_alignments = len(alignments)
    is_letter = np.core.defchararray.isalpha(alignments)
    frac_letters = is_letter.sum(0) * 1.0 / n_alignments
    return alignments[:, frac_letters >= frac_alpha_cutoff]
