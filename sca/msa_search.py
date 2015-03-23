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
    query = "".join([c for c in resname_unique if c.isalnum()])
    # Since non alphabetic chars in query are not used in the code
    scores = []

    print 'query: ' + query + ' of length ' + str(len(query))

    # Finding the sequence with the max score
    for index in range(len(alignments)):
        temp_align = alignments[index, np.core.defchararray.isalpha(alignments[index])]
        scores.append(smith_water_align(temp_align.tostring(), query).score)
    max_score = max(scores)
    max_score_index = scores.index(max_score)

    print 'alignments is of length ' + str(len(alignments))
    print 'scores: ' + str(scores)
    print 'max score: ' + str(max(scores))

    if truncate_or_not:
        is_letter = np.core.defchararray.isalpha(alignments[max_score_index])
        alignments_trunc = alignments[:, is_letter]
        sw_result = smith_water_align(alignments_trunc[max_score_index].tostring(), query)
        print '     truncated alignment using sequence #' + str(max_score_index) + \
            '  (score: ' + str(max_score) + ')'
    else:
        sw_result = smith_water_align(alignments[max_score_index].tostring(), query)
        print 'tophit is sequence #' + str(max_score_index) + '  (score: ' + str(max_score) + ')'

    print 'ref: ' + alignments[max_score_index].tostring()
    print 'query: ' + query
    top_score, best_align, query_start, ref_start, query_end, ref_end = sw_result.dump()

    for row in best_align:
        print row
    print 'length of best_align is: ' + str(len(best_align[0]))
    print top_score, query_start, ref_start, query_end, ref_end

    i = query_start - 1
    ats = [resseq_icode_unique[i]]
    for j in range(1, len(best_align[0])):
        # temp = best_align[1][j]
        # if temp == '|' or (temp == '.' and np.core.defchararray.isalpha(best_align[2][j])) or \
        #         (temp == ' ' and np.core.defchararray.isalpha(best_align[2][j])):
            # i += 1
            # ats.append(resseq_icode_unique[i])
        i += 1
        ats.append(resseq_icode_unique[i])
    assert len(best_align[0]) == len(best_align[2]), 'The output of smith_water_align i.e., ' \
                                                     'best_align has sequences of unequal lengths'
    profile1 = seq_profile(np.array([best_align[0], best_align[2]], 'c'))
    profile2 = seq_profile(alignments[max_score_index])
    prof, h1, h2 = profile_align(profile1, profile2)
    ats_new = [' '] * len(alignments[0])
    for i in range(len(h2)):
        if h2[i] in h1:
            ats_new[i] = ats[np.where(h1 == h2[i])[0]]
    if truncate_or_not:
        return max_score_index, ats_new, best_align, alignments_trunc
    else:
        return max_score_index, ats_new, best_align

def profile_align(prof1, prof2):
    """
    Profile-to-profile Needleman-Wunsch global alignment algorithm
    :param prof1: aminoacid-profile of a sequence (output of seqprofile function);
    must be a numpy array of 20 rows
    :param prof2: aminoacid-profile of a sequence (output of seqprofile function);
    must be a numpy array of 20 rows
    :return:
    """
    assert len(prof1) == 20 and len(prof2) == 20, 'The profiles prof1 and prof2 ' \
                                                  'must be numpy arrays with 20 rows'
    # The following are the assumptions made
    # extendgap = False
    # adjust_oldgap = True
    # adjust_endgap = False
    # v2_flag = False
    gap_penalty = -8
    gapgap_function = lambda a, b: 0.1 * a
    gapres_function = lambda a, b: 0.1 * b

    blosum50 = np.loadtxt('./sca/blosum50.txt')
    len1, len2 = len(prof1[0]), len(prof2[0])
    gap1, gap2 = np.ones(len1 + 1) * gap_penalty, np.ones(len2 + 1) * gap_penalty
    prof1 = np.concatenate((prof1, np.zeros((1, len(prof1[0])))))
    prof2 = np.concatenate((prof2, np.zeros((1, len(prof2[0])))))

    num_sym = 21
    scoring_matrix = np.zeros((num_sym, num_sym))
    scoring_matrix[0:num_sym - 1, 0:num_sym - 1] = blosum50[0:num_sym - 1, 0:num_sym - 1]
    sm = np.mean(np.diagonal(scoring_matrix))
    sx = sum(sum(scoring_matrix - np.diag(np.diagonal(scoring_matrix))))
    gapgap_const = gapgap_function(sm, sx)
    gapres_const = gapres_function(sm, sx)
    scoring_matrix[num_sym - 1, :] = gapres_const
    scoring_matrix[:, num_sym - 1] = gapres_const
    scoring_matrix[num_sym - 1, num_sym - 1] = gapgap_const

    gap_weight1, gap_weight2 = sum(prof1[:-1, :]), sum(prof2[:-1, :])
    temp1 = np.vstack((np.append(gap_weight1, max(gap_weight1)),
                       np.insert(gap_weight1, 0, max(gap_weight1))))
    gap1 = gap1 * np.min(temp1, axis=0)
    temp2 = np.vstack((np.append(gap_weight2, max(gap_weight2)),
                       np.insert(gap_weight2, 0, max(gap_weight2))))
    gap2 = gap2 * np.min(temp2, axis=0)

    f, pointer = needle_wunsch_align(prof1, prof2, scoring_matrix,
                                     gap1, gap2, gap_weight1, gap_weight2)

    i, j = len2 + 1, len1 + 1
    path = np.zeros((len1 + len2, 2))
    step = 1
    score = f[-1, -1]

    while i > 1 or j > 1:
        if pointer[i - 1, j - 1] == 1:
            i, j = i - 1, j - 1
            path[step - 1, :] = np.array([j, i])
        elif pointer[i - 1, j - 1] == 2:
            i -= 1
            path[step - 1, 1] = i
        elif pointer[i - 1, j - 1] == 4:
            j -= 1
            path[step - 1, 0] = j
        else:
            raise Exception
        step += 1
    path = path[:step - 1, :]
    path = np.flipud(path)
    prof = np.zeros((num_sym, step - 1))
    mask1 = path[:, 0] > 0
    mask2 = path[:, 1] > 0
    prof[:, mask1] = prof1
    prof[:, mask2] = prof[:, mask2] + prof2
    prof[num_sym - 1, ~mask1] = prof[num_sym - 1, ~mask1] + np.mean(sum(prof1))
    prof[num_sym - 1, ~mask2] = prof[num_sym - 1, ~mask2] + np.mean(sum(prof2))
    h1 = np.where(mask1 != 0)[0]
    h2 = np.where(mask2 != 0)[0]
    return prof, h1, h2


def needle_wunsch_align(prof1, prof2, scoring_matrix, gap1, gap2, gap_weight1, gap_weight2):
    """

    :param prof1:
    :param prof2:
    :param scoring_matrix:
    :param gap1:
    :param gap2:
    :param gap_weight1:
    :param gap_weight2:
    :return:
    """
    len1, len2 = len(prof1[0]), len(prof2[0])
    assert gap1.shape == (len1 + 1,), 'gap1 should be a numpy array of ' \
                                        '1 row, (len(prof1) + 1) cols'
    assert gap2.shape == (len2 + 1,), 'gap1 should be a numpy array of ' \
                                        '1 row, (len(prof2) + 1) cols'
    assert gap_weight1.shape == (len1,), 'gap_weight1 should be a numpy array ' \
                                               'of 1 row, len(prof1) cols'
    assert gap_weight2.shape == (len2,), 'gap_weight2 should be a numpy array ' \
                                               'of 1 row, len(prof2) cols'
    # f is the storage used (since it is a dynamic programming algo.)
    f = np.zeros((len2 + 1, len1 + 1))
    f[1:, 0] = gap1[0] * gap_weight2[0] * np.arange(1, len2 + 1)
    f[0, 1:] = gap2[0] * gap_weight1[0] * np.arange(1, len1 + 1)
    # pointer is the matrix used for backtracking
    pointer = np.ones((len2 + 1, len1 + 1)) * 4
    pointer[:, 0] = 2
    pointer[0, 0] = 1

    ptr = pointer[:, 1]
    curr_f_col = f[:, 0]
    scores = prof2.transpose(1, 0).dot(scoring_matrix.dot(prof1))

    for outer in np.arange(1, len1 + 1):
        last_f_col = curr_f_col
        curr_f_col = f[:, outer]
        best = curr_f_col[0]

        for inner in np.arange(1, len2 + 1):
            up = best + gap1[outer] * gap_weight2[inner - 1]
            left = last_f_col[inner] + gap2[inner] * gap_weight1[outer - 1]
            diagonal = last_f_col[inner - 1] + scores[inner - 1, outer - 1]

            if up > left:
                best = up
                pos = 2
            else:
                best = left
                pos = 4

            if diagonal >= best:
                best = diagonal
                ptr[inner] = 1
            else:
                ptr[inner] = pos
            curr_f_col[inner] = best

        f[:, outer] = curr_f_col
        pointer[:, outer] = ptr
    return f, pointer


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
