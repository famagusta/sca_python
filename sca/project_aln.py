'''This script contains functions to compute the projected
   weighted alignment from the alignment, the weighted 3d
   alignment tensor and the weight matrix'''

from helper import aa_freq, get_algn_shape, compute_pos_corr_matrix,\
    compute_seq_corr_matrix
import numpy as np
from numpy import linalg as LA
from numpy import zeros as zeros


class ProjectionMatrices(object):
    '''class that stores the claculated projection matrices
       and returns them'''
    def __init__(self, prj_wt_aln, proj_vect):
        self.prj_wt_aln = prj_wt_aln
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

    proj_vect = zeros((no_pos, no_aa))
    wt_mat = zeros((no_aa, no_pos))

    for i in range(0, no_pos):
        for aa in range(0, no_aa):
            wt_mat[aa, i] = weight[i, aa]*algn_aa_freq[aa, i]
        if LA.norm(wt_mat[:, i]) > 0:
            proj_vect[i, :] = wt_mat[:, i] / LA.norm(wt_mat[:, i])

    prj_wt_aln = zeros((no_seq, no_pos))

    for i in range(0, no_pos):
        for aa in range(0, no_aa):
            prj_wt_aln[:, i]\
                = prj_wt_aln[:, i]\
                + proj_vect[i, aa] * weighted_3d_algn[:, i, aa]
    return ProjectionMatrices(prj_wt_aln, proj_vect)
