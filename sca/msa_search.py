"""
Pairwise alignments between a query sequence read from a pdb file
and each sequence in a Multiple-Sequence Alignment
Main function msa_search with other helper functions
"""
import swalign
import numpy as np

AA = [['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
       'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
       'ASX', 'GLX', 'XAA', 'END', 'GAP'],
      ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
       'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*', '-']]
MATCH = 2
MISMATCH = -1
GAP_PENALTY = -8
GAP_EXTEND_PENALTY = -8


def msa_search(residue_list, alignments, truncate_or_not=False):
    """
    This function makes pairwise alignments between a query sequence (from
    the pdb file and chain ID) and every sequence in an MSA (alignment),
    finds the tophit sequence, and then attempts to make a residue number
    list (ats) that relates alignment numbering to structure numbering.  The
    function works with pdb files that comprise the standard pdb format
    according to the pdb.org.
    :param residue_list: Residue list obtained by reading a standard pdb file
    using read_pdb() function in sca_read.py file
    :param alignments: array consisting of alignments read from a .free or
    .fasta file
    :param truncate_or_not:
    :return:
    """
    # TODO: POSSIBLE BUG in MSAsearch.m (Octave version) in the first try-catch
    # block in line 72-73 indsort contains indices corresponding to atomnum_cat
    # indsort contains indices corresponding to atomnum_cat
    # but is used to select entries from atom_fields
    # atom_fields might need to be replaced by atomnum_fields in line 73
    # IF It's not a bug, make sure to change the following code.

    # icode is assumed (in sync with the get_structure method of PDBParser obj)
    # to be present as third entry in the tuple residue.id

    resseq_icode_list = \
        [(str(residue.id[1]) + residue.id[2]) for residue in residue_list]
    resname_list = [residue.resname for residue in residue_list]

    resseq_icode_unique, indices_unique = unique(resseq_icode_list)
    resname_unique = [resname_list[index] for index in indices_unique]

    # resseq_icode_unique ~ resnum_list and resname_unique ~ aa_model in the
    # Octave version

    # Replacing the 3-letter resnames in resname_unique by the
    # corresponding 1-letter resnames
    resname_unique =\
        [AA[1][AA[0].index(resname)] for resname in resname_unique]
    query = "".join(resname_list)
    scores = []
    for index in range(len(alignments)):
        scores.append(smith_water_align(query, alignments[index].tostring()))
    max_score = max(scores)
    max_score_index = scores.index(max_score)

    if truncate_or_not:
        is_letter = np.core.defchararray.isalpha(alignments[max_score_index])
        alignment_trunc = alignments[:, is_letter]
        print '     truncated alignment using sequence #' + \
            str(max_score_index) + \
            '  (score: ' + str(max_score) + ')'
    else:
        print 'tophit is sequence #' + str(max_score_index) + \
            '  (score: ' + str(max_score) + ')'

    if 'alignment_trunc' in locals():
        return max_score_index, alignment_trunc
    else:
        return max_score_index


def smith_water_align(seq1, seq2):
    """
    Applies Smith-Waterman sequence alignment algo. on the two input
    sequences and returns the score
    :param seq1: aminoacid sequence 1
    :param seq2: aminoacid sequence 2
    :param match: points for each match
    :param mismatch: points for each mismatch
    :param gap_penalty: points for each new gap
    :param gap_extension_penalty: points for each extension of gap
    :return:
    """
    scoring = swalign.NucleotideScoringMatrix(MATCH, MISMATCH)
    align_obj = swalign.LocalAlignment(scoring, GAP_PENALTY,
                                       GAP_EXTEND_PENALTY)
    align = align_obj.align(seq1, seq2)
    return align.score


def unique(my_list):
    """
    Returns a list made out of unique elements in my_list in the same
    order as they appear the first time in my_list
    :param my_list: input list with possible duplicate elements
    :return: a list with unique elements in the order they appear in the input
    list, a list with the indices of the unique elements in the input list
    """

    unique_elems = []
    indices = []
    for index in range(len(my_list)):
        item = my_list[index]
        if item not in unique_elems:
            unique_elems.append(item)
            indices.append(index)
    return unique_elems, indices
