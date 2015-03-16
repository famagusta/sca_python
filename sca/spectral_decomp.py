'''This script contains function definition for spectral decompostion.
    It also contains a class definition for returning eigenvalues &
    eigenvectors'''

from numpy import linalg as LA
from numpy import argsort as argsort
from numpy import sign as sign
from numpy import mean as mean
from numpy import zeros as zeros
from numpy import transpose as transpose
from numpy import dot as dot
from numpy import mean as mean
from numpy import random
from numpy import squeeze as squeeze
from helper import binrep, compute_pos_corr_matrix
from weighted_aln import weight_aln
from project_aln import project_aln
import numpy as np


class SpectralDecomp:
    '''This class stores eigen values & eigenvectors that result from the spectral
        decomposition of the sca matrix'''
    def __init__(self, C_rnd, pos_lbd, pos_ev, pos_lbd_rnd, pos_ev_rnd):
        self.C_rnd = C_rnd
        self.pos_lbd = pos_lbd
        self.pos_ev = pos_ev
        self.pos_lbd_rnd = pos_lbd_rnd
        self.pos_ev_rnd = pos_ev_rnd


class eig_pair:
    '''this class contains the eigenvalues and eigenvectors of a matrix'''
    def __init__(self, eig_vct, eig_val):
        self.eig_vct = eig_vct
        self.eig_val = eig_val


def spectral_decomp(sca_obj, no_samples, no_ev):
    '''This function computes the spectral (eigenvalue) decomposition of the
        correlation matrices in SCA. It carries out the N samples trial of such
        decomposition for randomized alignments.
        Input : SCA correlation matrices
        Output: Eigenvalues and Eigenvectors of actual & randomized alignments.
        QUESTION : Do we want to return randomized eignevectors?
    '''
    rnd_perm = random.permutation
    pwX = sca_obj.pwX
    shape = pwX.shape
    no_seq = shape[0]
    no_pos = shape[1]

    no_aa = 20

    # TODO : add a condition for no of arguments less that 2
    eig_obj = eig_vct(sca_obj.Cp, all_eig='all')
    ev = eig_obj.eig_vct
    lbd = eig_obj.eig_val

    # random trials to get randomly expected eigvalue spectra
    C_rnd = zeros((no_pos, no_pos))
    C_tmp = zeros((no_pos, no_pos, no_samples))

    algn_rnd = sca_obj.alignment
    print "Computing Random Trials . . . "
    for s in range(0, no_samples):
        for pos in range(0, no_pos):
            perm_seq = rnd_perm(no_seq)
            algn_rnd[:, pos] = sca_obj.alignment[perm_seq[:], pos]
        algn_rnd_3d_bin = binrep(algn_rnd)
        rnd_wt_matrices = weight_aln(algn_rnd_3d_bin)
        rnd_wt_3d_algn = rnd_wt_matrices.weighted_3d_algn

        prj_wt_aln_rnd = zeros((no_seq, no_pos))
        for i in range(0, no_pos):
            for a in range(0, no_aa):
                prj_wt_aln_rnd[:, i] = \
                    prj_wt_aln_rnd[:, i] + \
                    sca_obj.pm[i, a]*rnd_wt_3d_algn[:, i, a]

        C_tmp[:, :, s] = compute_pos_corr_matrix(prj_wt_aln_rnd, no_seq)
        C_rnd = C_rnd + squeeze(C_tmp[:, :, s])
        if s % 10 == 0:
            print "Iteration " + str(s)
    C_rnd = C_rnd/no_samples

    lbd_rnd = zeros((no_samples, no_pos))
    ev_rnd = zeros((no_samples, no_pos, no_ev))
    lbd_rnd_tmp = zeros((no_pos, 1))
    ev_rnd_tmp = zeros((no_pos, no_pos))

    for s in range(0, no_samples):
        eig_pair_temp = eig_vct(C_tmp[:, :, s], all_eig='all')
        lbd_rnd[s, :] = eig_pair_temp.eig_val
        ev_rnd[s, :, :] = eig_pair_temp.eig_vct[:, 0:no_ev]
    return SpectralDecomp(C_rnd, lbd, ev, lbd_rnd, ev_rnd)


def eig_vct(A, **options):
    ''' computes the first k (all if k>rank(A)) eigvectors of a matrix A
    and their corresponding eigvalues'''
    # v stores eigvalues d stores eigvectors in columns v[:,i]
    # corresponds to w[i]

    if options.get('all_eig') == 'all':
        v, d = LA.eig(A)
        return eig_pair(d, v)
    else:
        v, d = LA.eig(A)
        l = sorted(v, reverse=True)
        i = argsort(-v)
        w = v[:, i[1:k]]

        r = l[1:k]

        for n in range(0, k):
            w[:, n] = sign(mean(w[:, n]))*w[:, n]
        return eig_pair(r, w)
