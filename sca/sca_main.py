''' This script contains functions for performing sca for protein family
    given a multiple sequence alignment and a number of optional argtuments
'''

import numpy as np
from helper import binrep, get_algn_shape
from helper import compute_pos_corr_matrix, compute_seq_corr_matrix
from project_aln import project_aln
from weighted_aln import weight_aln
from numpy import matrix as matrix


class SCA:
    '''class structure that contains the output for the
        result of SCA calculations in the function before
    '''
    def __init__(self, alignment,
                 algn_3d_bin, weighted_3d_algn,
                 pwX, pm,
                 weight,  # para,
                 Cp, Cs):
        self.alignment = alignment
        self.algn_3d_bin = algn_3d_bin
        self.weighted_3d_algn = weighted_3d_algn
        self.pwX = pwX
        self.pm = pm
        self.weight = weight
        #  self.para = para
        self.Cp = Cp
        self.Cs = Cs


def sca(alignment):
    ''' main function that performs sca
        does not support user defined weighting function yet
        uses default LB entropy function and its derivative
        for weighting.
        output is a class containing
            alignment         : original input alignment
            algn_3d_bin       : binary alignment tensor (MxLx20)
            weighted_3d_algn  : MxLx20 weighted alignment tensor
            pwX               : the MxL projected weighted alignment matrix
            pm                : the Lx20 projection matrix
            weight            : the Lx20 weight matrix
            'weight: the weighting function handle - REMOVED only 1 weight fun
            para  : the value of the optional parameter
            Cp    : the SCA positional correlation matrix
            Cs    : the SCA sequence correlation matrix
        TODO : add these options in future if required by users
    '''

    algn_shape = get_algn_shape(alignment)
    no_aa = algn_shape.no_aa
    no_seq = algn_shape.no_seq
    no_pos = algn_shape.no_pos

    # Step 1 : Conversion of MSA into a 3D binary tensor X(a, i, s) with
    # x(a, i, s)=1 <=> sequence s has amino acid a at location i
    # call binrep to do this
    algn_3d_bin = binrep(alignment)

    # Step 2 : Calcualte weighted 3D tensor wX
    # call weight_aln to do this
    weight_matrices = weight_aln(algn_3d_bin)
    weight = weight_matrices.weight
    weighted_3d_algn = weight_matrices.weighted_3d_algn

    # Step 3 : Calculate 2D projection alignment matrix, pwX from 3D weighted
    # alignment tensor wX
    # call project_aln to do this
    proj_matrices = project_aln(alignment, weighted_3d_algn, weight)
    prj_wt_aln = matrix(proj_matrices.prj_wt_aln)
    proj_vect = matrix(proj_matrices.proj_vect)

    # Step 4 : Calculate SCA matrices

    pos_corr = compute_pos_corr_matrix(prj_wt_aln, no_seq)

    seq_corr = compute_seq_corr_matrix(prj_wt_aln, no_pos)

    return SCA(alignment, algn_3d_bin, weighted_3d_algn,
               prj_wt_aln, proj_vect, weight, pos_corr, seq_corr)
