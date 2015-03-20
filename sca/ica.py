"""
Independent component analysis (ICA) using maximum likelihood,
square mixing matrix and no noise (Infomax).
Source prior is assumed to be p(s)=1/pi*exp(-ln(cosh(s))).

Reference:
A. Bell and T.J. Sejnowski(1995).
An Information-Maximization Approach to Blind Separation
and Blind Deconvolution
Neural Computation, 7:1129-1159.

"""

from numpy import eye, fix, dot, exp, reshape, size, transpose


def ica(x, r, no_iter):
    '''computes basic infomac ICA
       inputs : x - LxM input matrix where L = #features
       and  M = #samples
       r : learning rate / relaxation parameter
       no_iter : number of iterations
    '''
    L = x.shape[0]
    M = x.shape[1]
    w = eye(L)  # initialization
    B = M  # not given the user to use lower B

    for i in range(0, no_iter-1):
        w_old = w
        j = 0
        for j in range(j, B, int(j+fix(M/B)*B)):
            u = dot(w, x[:, j:j+B])
            mid = dot(dot(B, eye(L)) + dot((1-2*(1/(1+exp(-u)))),
                                           transpose(u)), w)
            w = w + dot(r, mid)
        if(i % 1000 == 0):
            print 'iteration' + str(i)
    return w
