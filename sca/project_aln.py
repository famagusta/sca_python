'''This script contains functions to compute the projected
   weighted alignment from the alignment, the weighted 3d
   alignment tensor and the weight matrix'''

from helper import aa_freq, get_algn_shape
import numpy as np
from numpy import linalg as LA


class ProjectionMatrices(object):
    '''class that stores the claculated projection matrices
       and returns them'''
    def __init__(self, proj_wt_algn, proj_vect):
        self.proj_wt_algn = proj_wt_algn
        self.proj_vect = proj_vect


def project_aln(alignment, weighted_3d_algn, weight):
    '''computes projected weighted alignment based on the
       alignment in ascii format, the weighted 3D alginment
       tensor (weighted_3d_algn) & the weight matrix (weight)
       it does so by computing the frequencies first call (aa_freq)
    '''
    algn_shape = get_algn_shape(alignment)
    no_seq = algn_shape.no_seq
    no_pos = algn_shape.no_pos
    no_aa = algn_shape.no_aa

    algn_aa_freq = aa_freq(alignment)

    proj_vect = np.zeros((no_pos, no_aa))
    wght_mat = np.zeros((no_aa, no_pos))

    for i in range(0, no_pos):
        for aa in range(0, no_aa):
            wght_mat[aa, i] = weight[i, aa]*algn_aa_freq[aa, i]
        if LA.norm(wght_mat[:, i]) > 0:
            proj_vect[i, :] = wght_mat[:, i] / LA.norm(wght_mat[:, i])

    proj_wt_algn = np.zeros((no_seq, no_pos))

    for i in range(0, no_pos):
        for aa in range(0, no_aa):
            proj_wt_algn[:, i]\
                = proj_wt_algn[:, i]\
                + proj_vect[i, aa] * weighted_3d_algn[:, i, aa]
    return ProjectionMatrices(proj_wt_algn, proj_vect)
