"""Functions for reading input files"""
from Bio.PDB.PDBParser import PDBParser
from Bio import AlignIO
import numpy as np
from Bio import Alphabet


def read_free(filename):
<<<<<<< HEAD
    # TODO: there might be a library function in Biopython
    seqs = []
    for line in open(filename):
        (label, seq_str) = line.strip().split()
        seqs.append(list(seq_str))  # seq-str.ljust(SEQ_LENGTH)
=======
    """
    Reading .free files
    Input: filename - name of the file
    Output: ndarray
    TODO: there might be a library function in Biopython
    """
    seqs = []
    for line in open(filename):
        (_, seq_str) = line.strip().split()
        seqs.append(list(seq_str)) # seq-str.ljust(SEQ_LENGTH)
>>>>>>> 859c79122e7f6b273fec9321094f6a2a7c498cac
    return np.array(seqs, np.character)


def read_pdb(filename, chainid):
<<<<<<< HEAD
    # TODO: think of a better data structure for the output
=======
    """
    Reading .pdb file
    Input:
        filename - name of the file
        chainid
    Output: ndarray
    TODO: think of a better data structure for the output
    """
>>>>>>> 859c79122e7f6b273fec9321094f6a2a7c498cac
    structure = PDBParser().get_structure('PDB_STRUCT', filename)
    model = structure.child_list[0]
    chain = model.child_dict[chainid]
    residue_list = chain.child_list
    return residue_list


def read_fasta(filename):
<<<<<<< HEAD
    msa = AlignIO.read(filename, 'fasta',
                       alphabet=Alphabet.Gapped(Alphabet.IUPAC.protein))
=======
    """
    Reading .fasta files
    Input: filename - name of the file
    Output: ndarray
    """
    msa = AlignIO.read(filename, 'fasta', alphabet = Alphabet.Gapped(Alphabet.IUPAC.protein))
>>>>>>> 859c79122e7f6b273fec9321094f6a2a7c498cac
    return np.array([list(rec) for rec in msa], np.character)
