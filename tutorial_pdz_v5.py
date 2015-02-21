'''
Analysis of an alignment of the PDZ domain family
Calculate SVD of a matrix using numpy functions
e.g.
A = floor(random.rand(4,4)*20-10)
linalg.svd(A)
for sparse matrices, there are methods in SciPy
'''

from sca import sca_read, msa_search, tutorial_helpers as helpers

PDB_FILE = './Inputs/1BE9.pdb'
FREE_FILE = './Inputs/al_pdz.free'
CHAINID = 'A'
FRAC_ALPHA_CUTOFF = 0.8


def tutorial():
    '''
    Analysis of an alignment of the PDZ domain family
    '''
    pdb = sca_read.read_pdb(PDB_FILE, CHAINID)
    alignments = sca_read.read_free(FREE_FILE)
    # n_alignments = len(alignments)
    # len_alignment = len(alignments[0])

    # truncate alignments to sequence positions with
    # gap frequency no greater than 20% - to avoid over-representation of gaps
    alignments = helpers.truncate(alignments, FRAC_ALPHA_CUTOFF)

    return msa_search.msa_search(pdb, CHAINID, alignments)


OUTPUT = tutorial()
