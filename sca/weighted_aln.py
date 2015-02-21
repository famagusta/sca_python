''' This script calculates the weighted alignment tensor and
the weight matrix using the binary tensor
representation of a multiple sequence alignment.
The default weighing function is the gradient of the relative
entropy function. Here we use the Kulback-Leibler relative entropy
function only'''

import numpy as np
import numpy.matlib
import math


class WeightMatrices(object):
    def __init__(self, weight, weighted_3d_algn):
        self.weight = weight
        self.weighted_3d_algn = weighted_3d_algn


def weight_aln(algn_3d_bin):
    '''The main function to calculate the weighted
    alignment tensor from the binary tensor representation
    of the multiple sequence alignment.'''

    # Constants ... the background frequencies for amino acids
    aa_freq = [0.073, 0.025, 0.050, 0.061, 0.042, 0.072, 0.023, 0.053,
               0.064, 0.089, 0.023, 0.043, 0.052, 0.040, 0.052, 0.073,
               0.056, 0.063, 0.013, 0.033]

    # MAIN Program
    shape = algn_3d_bin.shape
    no_seq = shape[0]
    no_pos = shape[1]
    no_aa = shape[2]

    # Convert the 3D binary tensor version of the MSA into a matrix
    # representation for computational simplicity
    algn_3d_mat = np.zeros((no_seq, no_aa*no_pos))
    algn_3d_mat = np.reshape(np.transpose(algn_3d_bin,
                             (2, 1, 0)), (no_aa*no_pos, no_seq), order='F')

    # Compute the frequencies and compute the weights given the reference
    # frequencies of amino acids, a weight function, and an optional parameter
    freq_vect = (np.transpose(algn_3d_mat).sum(0))/no_seq
    weight_vect = np.ones((1, no_aa*no_pos))
    aa_rep_vect = numpy.matlib.repmat(aa_freq, 1, no_pos)

    for i in range(0, no_aa*no_pos):
        weight_vect[0, i] = deriv_entropy(freq_vect[i],
                                          aa_rep_vect[0, i])

    weight = np.zeros((no_pos, no_aa))

    for i in range(0, no_pos):
        for aa in range(0, no_aa):
            weight[i, aa] = weight_vect[0, no_aa*i + aa]

    weighted_3d_algn = np.zeros((no_seq, no_pos*no_aa))
    weight_reshaped = np.reshape(weight, (1, no_pos, no_aa), order='F')
    weighted_3d_algn = np.tile(weight_reshaped,
                               (no_seq, 1, 1))*algn_3d_bin

    return WeightMatrices(weight, weighted_3d_algn)


def deriv_entropy(f, q):
    '''Derivative of relative entropy dD(f||q)/df'''

    D = 0
    if f > 0 and f < 1:
        D = abs(math.log(f*(1-q)/(q*(1-f))))
    return D
