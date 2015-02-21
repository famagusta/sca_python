import sca_read
import swalign

AA = [['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ASX', 
    'GLX', 'XAA', 'END', 'GAP'],
    ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 
    'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*', '-']]
MATCH = 2
MISMATCH = -1
GAP_PENALTY = -8
GAP_EXTEND_PENALTY = -8

def msa_search(residue_list, chainid, alignments, truncate_or_not = False):
    ''''''
    # TODO: POSSIBLE BUG in MSAsearch.m in the first try-catch block in line 72-73
    # indsort contains indices corresponding to atomnum_cat but is used to select entries from atom_fields
    # atom_fields might need to be replaced by atomnum_fields in line 73
    # IF It's not a bug, make sure to change the following code.
    
    # icode is assumed (in sync with the get_structure method of PDBParser obj.) 
    # to be present as third entry in the tuple residue.id
    resseq_icode_list = [(str(residue.id[1]) + residue.id[2]) for residue in residue_list]
    resname_list = [residue.resname for residue in residue_list]
    
    resseq_icode_unique, indices_unique = unique(resseq_icode_list)
    resname_unique = [resname_list[index] for index in indices_unique]
    # resseq_icode_unique ~ resnum_list and resname_unique ~ aa_model in the Octave version
    
    # Replacing the 3-letter resnames in resname_unique by the corresponding 1-letter resnames
    resname_unique = [AA[1][AA[0].index(resname)] for resname in resname_unique]
    query = "".join(resname_list)
    scores = []
    for index in range(len(alignments)):
        scores.append(smith_water_align(query, alignments[index], MATCH, MISMATCH, GAP_PENALTY, GAP_EXTEND_PENALTY)
    max_score = max(scores)

def smith_water_align(seq1, seq2, match, mismatch, gap_penalty, gap_entension_penalty):
    scoring = swalign.IdentityScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty, gap_extension_penalty)
    align = sw.align(seq1, seq2)
    return align.score

def unique(my_list):
    '''
    Returns a list made out of unique elements in my_list in the same
    order as they appear the first time in my_list
    '''
    unique_elems = []
    indices = []
    for index in range(len(my_list)):
        item = my_list[index]
        if item not in unique_elems:
            unique_elems.append(item)
            indices.append(index)
    return unique_elems, indices
