''' Convert standard ASCII MSA MxL (M~sequences, L~positions) to
    a MxLx20 3D binary tensor (X3d)'''
import numpy as np


def binrep(algn):
    '''Main function that creates the 3D binary tensor'''
    shape_algn = algn[1].shape
    N_aa = 20
    alg = lett2num(algn)
    shape = [shape_algn[0], shape_algn[1], N_aa]
    X3d = np.zeros(shape)
    for a in range(0, N_aa):
        X3d[:, :, a] = (alg == a)
    return X3d


def lett2num(msa_lett):
    '''translates an alignment from a character representation of
    amino acids ACDEFGHIKLMNPQRSTVWY to a numeric representation [1-20]'''
    code = 'ACDEFGHIKLMNPQRSTVWY'

    shape = msa_lett[1].shape
    msa_num = np.zeros(shape)

    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            lett = msa_lett[1][i][j]
            for index in range(0, len(code)):
                if lett == code[index]:
                    msa_num[i, j] = index

    return msa_num
