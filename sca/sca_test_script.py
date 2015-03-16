from sca import *
from sca_read import *
from helper import *
from weighted_aln import *
from project_aln import *
from numpy import linalg as LA
from spectral_decomp import *
import numpy as np


algn = read_free(
    '/Users/robinphilip/Documents/TAPAS/tapas_app/Inputs/al_pdz.free')
sca_algn = sca(algn)
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
