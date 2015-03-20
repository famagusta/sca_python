from sca import *
from sca_read import *
from helper import *
from weighted_aln import *
from project_aln import *
from numpy import linalg as LA
from scipy.stats import lognorm
from scipy.stats import t
from spectral_decomp import *
import numpy as np
from numpy import arange
from numpy import log
from sklearn.decomposition import FastICA


base_input_dir = '/Users/robinphilip/Documents/TAPAS/tapas_app/Inputs/'

pdz_file = base_input_dir + 'al_pdz.free'
sprot_file = base_input_dir + 'al_S1A_1388.free'

uncorr_algo = 'svd'

distbn_to_fit = 'lognormal'

algn = read_free(sprot_file)
sca_algn = sca(algn)
algn_shape = get_algn_shape(algn)
no_pos = algn_shape.no_pos
no_seq = algn_shape.no_seq
no_aa = algn_shape.no_aa
print 'Testing SCA module :'
print 'algn_3d_bin hash :' + str(np.sum(np.square(sca_algn.algn_3d_bin)))
print 'weighted_3d_algn hash :' +\
    str(np.sum(np.square(sca_algn.weighted_3d_algn)))
print 'weight hash : ' + str(np.sum(np.square(sca_algn.weight)))
print 'pwX hash : ' + str(np.sum(np.square(sca_algn.pwX)))
print 'pm hash : ' + str(np.sum(np.square(sca_algn.pm)))
print 'Cp has : ' + str(np.sum(np.square(sca_algn.Cp)))
print 'Cs hash : ' + str(np.sum(np.square(sca_algn.Cs)))

spect = spectral_decomp(sca_algn, 100, 100)
print 'spect lb hash : ' + str(np.sum(np.square(spect.pos_lbd)))
print 'spect ev hash : ' + str(np.sum(np.square(spect.pos_ev)))
print 'spect ldb_rnd hash : ' + str(np.sum(np.square(spect.pos_lbd_rnd)))
print 'spect ev hash : ' + str(np.sum(np.square(spect.pos_ev_rnd)))

svd_output = LA.svd(sca_algn.pwX)
U = svd_output[0]
sv = svd_output[1]
V = svd_output[2]

# calculate the matrix Pi = U*V'
# this provides a mathematical mapping between
# positional and sequence correlation

n_min = min(no_seq, no_pos)
Pi = dot(U[:, 0:n_min-1], transpose(V[:, 0:n_min-1]))
U_p = dot(Pi, spect.pos_ev)

pd = t.fit(spect.pos_ev[:, 0], floc=0)
# floc = 0 holds location to 0 for fitting
print pd

ica = FastICA(max_iter=1000000, algorithm='parallel', tol=0.0001)


def test_z(filename, uncorr_algo, distbn_to_fit):
    '''test case for pdz domain proteins'''
    algn = read_free(filename)
    sca_algn = sca(algn)
    algn_shape = get_algn_shape(algn)
    no_pos = algn_shape.no_pos
    no_seq = algn_shape.no_seq
    no_aa = algn_shape.no_aa
    print 'Testing SCA module :'
    print 'algn_3d_bin hash :' + str(np.sum(np.square(sca_algn.algn_3d_bin)))
    print 'weighted_3d_algn hash :' +\
        str(np.sum(np.square(sca_algn.weighted_3d_algn)))
    print 'weight hash : ' + str(np.sum(np.square(sca_algn.weight)))
    print 'pwX hash : ' + str(np.sum(np.square(sca_algn.pwX)))
    print 'pm hash : ' + str(np.sum(np.square(sca_algn.pm)))
    print 'Cp has : ' + str(np.sum(np.square(sca_algn.Cp)))
    print 'Cs hash : ' + str(np.sum(np.square(sca_algn.Cs)))

    spect = spectral_decomp(sca_algn, 100, 100)
    print 'spect lb hash : ' + str(np.sum(np.square(spect.pos_lbd)))
    print 'spect ev hash : ' + str(np.sum(np.square(spect.pos_ev)))
    print 'spect ldb_rnd hash : ' + str(np.sum(np.square(spect.pos_lbd_rnd)))
    print 'spect ev hash : ' + str(np.sum(np.square(spect.pos_ev_rnd)))

    svd_output = LA.svd(sca_algn.pwX)
    U = svd_output[0]
    sv = svd_output[1]
    V = svd_output[2]

    # calculate the matrix Pi = U*V'
    # this provides a mathematical mapping between
    # positional and sequence correlation

    n_min = min(no_seq, no_pos)
    Pi = dot(U[:, 0:n_min-1], transpose(V[:, 0:n_min-1]))
    U_p = dot(Pi, spect.pos_ev)

    distbn = get_distbn(distbn_to_fit)
    pd = distbn.fit(spect.pos_ev[:, 0], floc=0)
    # floc = 0 holds location to 0 for fitting
    print pd

    p_cutoff = 0.8  # cutoff for the cdf
    xhist = arange(0, 0.4, 0.01)
    x_dist = arange(min(xhist), max(xhist), (max(xhist) - min(xhist))/100)
    cdf = lognorm.cdf(x_dist, pd[0], pd[1], pd[2])
    # Use case : lognorm.cdf(x, shape, loc, scale)

    jnk = min(abs(cdf - p_cutoff))
    x_dist_pos_right = np.argmin(abs(cdf-p_cutoff))
    cutoff_ev = x_dist[x_dist_pos_right]
    sector_def = np.array(np.where(spect.pos_ev[:, 0] > cutoff_ev)[0])[0]


def get_distbn(distbn_to_fit):
    '''chooses the distribution from scipy based on user input'''
    if(distbn_to_fit == 'tlocscale'):
        distbn = t
    else:
        distbn = lognorm
