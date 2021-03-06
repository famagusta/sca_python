'''Contains a bunch of functions to make life easier'''

import numpy as np
from numpy import transpose as tp
from numpy import mean as mean
from numpy import dot as dot
from numpy import abs as npabs
from numpy import matrix as matrix


class AlignmentShape(object):
    '''this class contains dimensions of alignment'''
    # assuming only 20 amino acids
    no_aa = 20

    def __init__(self, no_pos, no_seq):
        self.no_pos = no_pos
        self.no_seq = no_seq


def get_algn_shape(algn):
    '''returns the dimensions of the algn
       it is a pain to extract this again and again'''
    shape = algn.shape
    return AlignmentShape(shape[1], shape[0])


def binrep(alignment):
    '''Main function that creates the 3D binary tensor
    Convert standard ASCII MSA MxL (M~sequences, L~positions) to
    a MxLx20 3D binary tensor (X3d)'''
    shape_algn = alignment.shape
    no_aa = 20
    num_alignment = lett2num(alignment)
    shape = [shape_algn[0], shape_algn[1], no_aa]
    algn_3d_bin = np.zeros(shape)
    for index in range(0, no_aa):
        algn_3d_bin[:, :, index] = (num_alignment == index + 1)
    return algn_3d_bin


def lett2num(alignment):
    '''translates an alignment from a character representation of
    amino acids ACDEFGHIKLMNPQRSTVWY to a numeric representation [1-20]
    0 is for representing gaps'''
    code = 'ACDEFGHIKLMNPQRSTVWY'

    shape = alignment.shape
    num_alignment = np.zeros(shape)

    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            lett = alignment[i][j]
            for index in range(1, len(code) + 1):
                if lett == code[index - 1]:
                    num_alignment[i, j] = index

    return num_alignment


def aa_freq(alignment):
    '''compute the frequency of each amino acid in a MSA
       at a given location
       Input  : alignment of dim MxLx20
       Output : matrix of size 20xL
    '''
    shape_obj = get_algn_shape(alignment)
    no_aa = shape_obj.no_aa
    no_pos = shape_obj.no_pos
    no_seq = shape_obj.no_seq

    algn_3d_bin = binrep(alignment)

    freq = np.zeros((no_aa, no_pos))

    for aa in range(0, no_aa):
        freq[aa, :] = algn_3d_bin[:, :, aa].sum(0)/no_seq

    return freq


def compute_pos_corr_matrix(X, N):
    '''computes positional correlation matrix'''
    pos_mat_prod = dot(tp(X), X)/N
    pos_avg_prod = dot(tp(matrix(mean(X, 0))), matrix(mean(X, 0)))
    pos_corr = npabs(pos_mat_prod - pos_avg_prod)
    return pos_corr


def compute_seq_corr_matrix(X, N):
    '''computes sequnece correlation matrix'''
    seq_mat_prod = dot(X, tp(X))/N
    seq_avg_prod = \
        dot(tp(matrix(mean(tp(X), 0))), matrix(mean(tp(X), 0)))

    seq_corr = npabs(seq_mat_prod - seq_avg_prod)
    return seq_corr
