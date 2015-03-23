from sca import *
from sca_read import *
from helper import *
from weighted_aln import *
from project_aln import *
from numpy import linalg as LA
from scipy.stats import lognorm, t
from spectral_decomp import *
import numpy as np
from numpy import arange, log, square, zeros, percentile, ptp,\
    histogram, argmax, argmin, array, where
from numpy import sum as mat_sum
from ica import *


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


# perform independent components calculations
kmax = 8
learnrate = 0.0001
iterations = 20000
w = ica(transpose(spect.pos_ev[:, 0:kmax]), learnrate, iterations)
ic_P = transpose(dot(w, transpose(spect.pos_ev[:, 0:kmax])))

print "ic_P hash :" + str(mat_sum(square(ic_P)))
# calculate the matrix Pi = U*V'
# this provides a mathematical mapping between
# positional and sequence correlation

n_min = min(no_seq, no_pos)
Pi = dot(U[:, 0:n_min-1], transpose(V[:, 0:n_min-1]))
U_p = dot(Pi, spect.pos_ev)


p_cutoff = 0.9
nfit = 3
cutoffs = zeros((nfit, 1))
sector_def = []

for i in range(0, nfit):
    nu, mu, sigma = t.fit(ic_P[:, i])
    q75, q25 = percentile(ic_P[:, i], [75, 25])
    iqr = q75 - q25
    binwidth = 2*iqr*pow(size(ic_P[:, i]), -1/3.0)  # Freedman-Diaconis rule
    nbins = round(ptp(ic_P[:, i])/binwidth)
    yhist, xhist = histogram(ic_P[:, i], nbins)
    x_dist = arange(min(xhist), max(xhist), (max(xhist) - min(xhist))/100)
    cdf_jnk = t.cdf(x_dist, nu, mu, sigma)
    pdf_jnk = t.pdf(x_dist, nu, mu, sigma)
    maxpos = argmax(pdf_jnk)
    tail = zeros((1, size(pdf_jnk)))
    if abs(max(ic_P[:, i])) > abs(min(ic_P[:, i])):
        tail[:, maxpos:] = cdf_jnk[maxpos:]
    else:
        tail[0:maxpos] = cdf_jnk[0:maxpos]
    x_dist_pos = argmin(abs(tail - p_cutoff))
    cutoffs[i] = x_dist[x_dist_pos]
    sector_def.append(array(where(ic_P[:, i] > cutoffs[i])[0])[0])
print sector_def


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
