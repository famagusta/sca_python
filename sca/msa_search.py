"""
Pairwise alignments between a query sequence read from a pdb file
and each sequence in a Multiple-Sequence Alignment
Main function msa_search with other helper functions
"""
import swalign
import numpy as np

AA = [['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
       'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ASX',
       'GLX', 'XAA', 'END', 'GAP'],
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
    function works with pdb files that comprise the standard pdb format according to the pdb.org.
    :param residue_list: Residue list obtained by reading a standard pdb file
    using read_pdb() function in sca_read.py file
    :param alignments: array consisting of alignments read from a .free or .fasta file
    :param truncate_or_not: boolean variable indicating whether to truncate alignment or not
    :return:
    """

    # icode is assumed (in sync with the get_structure method of PDBParser obj.)
    # to be present as third entry in the tuple residue.id
    resseq_icode_list = [(str(residue.id[1]) + residue.id[2]) for residue in residue_list]
    resname_list = [residue.resname for residue in residue_list]

    # Removing duplicates
    resseq_icode_unique, indices_unique = unique(resseq_icode_list)
    resname_unique = [resname_list[index] for index in indices_unique]
    # resseq_icode_unique ~ resnum_list and resname_unique ~ aa_model in the Octave version

    # Replacing the 3-letter resnames in resname_unique by the corresponding 1-letter resnames
    resname_unique = [AA[1][AA[0].index(resname)] for resname in resname_unique]
    print resname_unique
    ref = "".join([c for c in resname_unique if c.isalnum()])
    # Since non alphabetic chars in ref are not used in the code
    scores = []

    print 'ref: ' + ref + ' of length ' + str(len(ref))

    # Finding the sequence with the max score
    for index in range(len(alignments)):
        temp_align = alignments[index, np.core.defchararray.isalpha(alignments[index])]
        scores.append(smith_water_align(temp_align.tostring(), ref).score)
    max_score = max(scores)
    max_score_index = scores.index(max_score)

    print 'alignments is of length ' + str(len(alignments))
    print 'scores: ' + str(scores)
    print 'max score: ' + str(max(scores))

    if truncate_or_not:
        is_letter = np.core.defchararray.isalpha(alignments[max_score_index])
        alignments_trunc = alignments[:, is_letter]
        sw_result = smith_water_align(alignments_trunc[max_score_index].tostring(), ref)
        print '     truncated alignment using sequence #' + str(max_score_index) + \
            '  (score: ' + str(max_score) + ')'
    else:
        sw_result = smith_water_align(alignments[max_score_index].tostring(), ref)
        print 'tophit is sequence #' + str(max_score_index) + '  (score: ' + str(max_score) + ')'

    print 'query: ' + alignments[max_score_index].tostring()
    print 'ref: ' + ref
    top_score, best_align, query_start, ref_start, query_end, ref_end = sw_result.dump()

    for row in best_align:
        print row
    print 'length of best_align is: ' + str(len(best_align[0]))
    print top_score, query_start, ref_start, query_end, ref_end

    i = ref_start
    ats = [resseq_icode_unique[i]]
    for j in range(1, len(best_align[0])):
        temp = best_align[1][j]
        if temp == '|' or temp == '.' or \
                (temp == ' ' and np.core.defchararray.isalpha(best_align[2][j])):
            i += 1
            ats.append(resseq_icode_unique[i])
    assert len(best_align[0]) == len(best_align[2]), 'The output of smith_water_align i.e., ' \
                                                     'best_align has sequences of unequal lengths'
    profile1 = seq_profile(np.array([best_align[0], best_align[2]], 'c'))
    profile2 = seq_profile(alignments[max_score_index])
    return profile1, profile2


def seq_profile(seqs, alphabet='aa', gaps='none', counts=False,
                count_ambiguous=False, limits=(0, float('inf'))):
    """
    Computes the sequence profile of a multiple alignment
    :param seqs: multiple sequence alignment with a sequence in each row - MUST be a numpy array
    :return: profile: array of size 20 (or 4 for Nucleotide sequences) X number of sequences
    in seqs, each row corresponding to an amino acid, each column to a position in the
    sequence, with the frequency of occurrence of amino acids in seqs as entries
    """
    assert alphabet in ['aa', 'nt', 'none'], 'Invalid input: alphabet'
    assert counts in [True, False], 'Invalid input: counts'
    assert gaps in ['all', 'noflanks', 'none'], 'Invalid input: gaps'
    assert count_ambiguous in [True, False], 'Invalid input: ambiguous'
    assert len(limits) == 2 and limits[0] < limits[1], 'Invalid input: limits'

    seqs = np.char.upper(seqs)
    seqs[seqs == '.'] = '-'
    seqs[seqs == ' '] = '-'
    # seqs must be a 2D array for our computations
    if len(seqs.shape) == 1:
        seqs = np.array([list(seqs)])
    num_seqs, num_cols = len(seqs), len(seqs[0])

    if gaps == 'noflanks':
        for seq_index in range(len(seqs)):
            gap_indices = np.nonzero(seqs[seq_index] != '-')
            seqs[seq_index, 0:gap_indices[0]] = ' '
            seqs[seq_index, gap_indices[len(gap_indices) - 1] + 1:] = ' '

    limits = max(0, limits[0]), min(limits[1], num_cols - 1)
    if limits[0] != 0 or limits[1] != num_cols - 1:
        seqs = seqs[:, limits[0]:limits[1]+1]
        num_cols = len(seqs[0])

    if alphabet == 'aa':
        # len(AA[0]) gives the number of standard AAs; we need to account for unknown AAs also
        occurrences = np.zeros((len(AA[0]) + 1, len(seqs[0])))
        seqs_in_ints = aa2int(seqs, unknown=25)
        for col in range(num_cols):
            occurrences[:, col] = np.bincount(seqs_in_ints[:, col],
                                              weights=np.ones(num_seqs), minlength=len(AA[0]) + 1)
        if gaps == 'none':
            freq_occur = occurrences[0:20, :]
        else:
            freq_occur = occurrences[np.concatenate((np.arange(0, 20), np.array([24]))), :]
        if count_ambiguous:
            freq_occur[0:20, :] = freq_occur[0:20, :] + \
                                  np.tile(occurrences[22, :]/20.0, np.array([20, 1]))
            freq_occur[5:7, :] = freq_occur[5:7, :] + \
                                 np.tile(occurrences[21, :]/2.0, np.array([2, 1]))
            freq_occur[2:4, :] = freq_occur[2:4, :] + \
                                 np.tile(occurrences[20, :]/2.0, np.array([2, 1]))
    elif alphabet == 'nt':
        raise Exception('seqprofile is not implemented as of now for Nucleotide sequences')
    else:
        raise Exception('seqprofile is not implemented as of now for sequences '
                        'other than Amino Acid sequences')

    if not counts:
        freq_sum = np.sum(freq_occur, axis=0)
        freq_sum[freq_sum == 0] = 1
        freq_occur = freq_occur / np.tile(freq_sum, np.array([len(freq_occur), 1]))

    return freq_occur


def aa2int(aas, unknown=0):
    """
    This function converts a list of strings of amino acids from letters to numbers
    :param aas: must be an array of characters in which each row corresponds
    to an amino-acid sequence
    :param unknown: the number used to represent unknown amino acids
    :return: numpy array of integers
    """
    aas = np.char.lower(aas)

    if len(aas) == 0:
        return np.array([])

    ints = np.array([1, 21, 5, 4, 7, 14, 8, 9, 10, unknown, 12, 11, 13, 3, unknown,
                     15, 6, 2, 16, 17, unknown, 20, 18, 23, 19, 22, 24, 25, unknown])
    # list_of_aas = np.array(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
    # 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '*', '-', '?'])

    aas = np.vectorize(ord)(aas) - 97  # convert characters into integers using ASCII code
    gaps = (aas == -52)
    stops = (aas == -55)
    aas[aas < 0] = 28
    aas[aas > 25] = 28
    aas[gaps] = 27
    aas[stops] = 26

    return ints[aas]


def smith_water_align(seq1, seq2):
    """
    Applies Smith-Waterman sequence alignment algo. on the two input sequences and returns the score
    :param seq1: aminoacid sequence 1
    :param seq2: aminoacid sequence 2
    :return:
    """
    scoring = swalign.NucleotideScoringMatrix(MATCH, MISMATCH)
    align_obj = swalign.LocalAlignment(scoring, GAP_PENALTY, GAP_EXTEND_PENALTY)
    align = align_obj.align(seq1, seq2)
    return align


def unique(my_list):
    """
    Returns a list made out of unique elements in my_list in the same
    order as they appear the first time in my_list
    :param my_list: input list with possible duplicate elements
    :return: a list with unique elements in the order they appear in the input list,
     a list with the indices of the unique elements in the input list
    """

    unique_elems = []
    indices = []
    for index in range(len(my_list)):
        item = my_list[index]
        if item not in unique_elems:
            unique_elems.append(item)
            indices.append(index)
    return unique_elems, indices
