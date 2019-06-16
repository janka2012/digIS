import logging
import os

from definitions import ROOT_DIR

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


if __name__ == "__main__":

    genbank_file = os.path.join(ROOT_DIR, "test_data", "GCA_000006805.1_ASM680v1_genomic.gbff")
    gb_list = read_gb(genbank_file)
    from pprint import pprint
    pprint(len(gb_list))
