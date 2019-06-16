import logging
import os


from Bio import SeqIO


def read_gb(filename):
    genbank_dict = {}
    if os.path.exists(filename):
        for genome in SeqIO.parse(filename, "genbank"):
            genbank_dict[genome.id] = genome.features
    else:
        logging.error("Filename {} does not exist.".format(filename))
        raise FileNotFoundError("{} file does not exist.".format(filename))

    return genbank_dict

