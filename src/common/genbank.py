import logging
import os

from Bio import SeqIO


def read_gb(filename):
    out_list = []
    if os.path.exists(filename):
        for rec in SeqIO.parse(filename, "genbank"):
            for f in rec.features:
                out_list.append(f)
    else:
        logging.error("Filename {} does not exist.".format(filename))
        raise FileNotFoundError("{} file does not exist.".format(filename))

    return out_list
