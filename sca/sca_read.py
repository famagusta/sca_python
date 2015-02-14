from Bio.PDB.PDBParser import PDBParser
import numpy as np

LABEL_LENGTH = 8
SEQ_LENGTH = 129

def read_free(filename):
    labels = []
    seqs = []
    for line in open(filename):
        (label, seq) = line.strip().split()
        labels.append(list(label.ljust(LABEL_LENGTH)))
        seqs.append(list(seq.ljust(SEQ_LENGTH)))
    return (np.array(labels), np.array(seqs))

def read_pdb(filename):
    #TODO: this is just a draft
    parser = PDBParser(filename)
    strcture = parser.get_structure()
    return structure

def read_fasta(filename):
    return [0]
