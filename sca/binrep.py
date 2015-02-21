''' Convert standard ASCII MSA MxL (M~sequences, L~positions) to
    a MxLx20 3D binary tensor (X3d)'''
import numpy as np


def binrep(alignment):
    '''Main function that creates the 3D binary tensor'''
    shape_algn = alignment[1].shape
    no_aa = 20
    num_alignment = lett2num(alignment)
    shape = [shape_algn[0], shape_algn[1], no_aa]
    algn_3d_bin = np.zeros(shape)
    for index in range(0, no_aa):
        algn_3d_bin[:, :, index] = (num_alignment == index)
    return algn_3d_bin


def lett2num(alignment):
    '''translates an alignment from a character representation of
    amino acids ACDEFGHIKLMNPQRSTVWY to a numeric representation [1-20]
    0 is for representing gaps'''
    code = 'ACDEFGHIKLMNPQRSTVWY'

    shape = alignment[1].shape
    num_alignment = np.zeros(shape)

    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            lett = alignment[1][i][j]
            for index in range(0, len(code)):
                if lett == code[index]:
                    num_alignment[i, j] = index

    return num_alignment
