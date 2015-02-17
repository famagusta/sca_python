''' This script contains functions for performing sca for protein family
    given a multiple sequence alignment and a number of optional argtuments
'''

import numpy as np
import binrep


class sca:
    '''class structure that contains the output for the
        result of SCA calculations in the function before
    '''
    def __init__(self, X3d, wX, pwX, pm, W, weight, para, Cp, Cs):
        self.X3d = X3d
        self.wX = wX
        self.pwX = pwX
        self.pm = pm
        self.W = W
        self.weight = weight
        self.para = para
        self.Cp = Cp
        self.Cs = Cs


def sca(algn):
    ''' main function that performs sca
        does not support user defined weighting function yet
        uses default LB entropy function and its derivative
        for weighting.
        output is a class containing
            X3d   : binary alignment tensor (MxLx20)
            wX    : MxLx20 weighted alignment tensor
            pwX   : the MxL projected weighted alignment matrix
            pm    : the Lx20 projection matrix
            W     : the Lx20 weight matrix
            weight: the weighting function handle
            para  : the value of the optional parameter
            Cp    : the SCA positional correlation matrix
            Cs    : the SCA sequence correlation matrix
        TODO : add these options in future if required by users
    '''
    # Step 1 : Conversion of MSA into a 3D binary tensor X(a, i, s) with
    # x(a, i, s)=1 <=> sequence s has amino acid a at location i
    # call binrep to do this
    X3d = binreq(algn)
    
    # Step 2 : Calcualte weighted 3D tensor wX
    # call weight_aln to do this

    # Step 3 : Calculate 2D projection alignment matrix, pwX from 3D weighted
    # alignment tensor wX
    # call project_aln to do this

    # Step 4 : Calculate SCA matrices
    return 0
