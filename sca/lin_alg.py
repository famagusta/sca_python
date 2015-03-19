'''this script contains functions for computing matrix decompositions
   - SVD : Singular Value Decomposition
   - ICA : Independent Component Analysis
'''

from numpy import linalg as LA


class SVD:
    '''this class holds the results of singular value decomposition
       this makes it easier to pass around in the code
    '''
    def __init__(self, u, s, v):
        self.u = u
        self.s = s
        self.v = v


def compute_svd(A):
    return 0
