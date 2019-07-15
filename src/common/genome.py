import logging
import os

from collections import OrderedDict
from Bio import SeqIO


class Genome:

    def __init__(self, genome_rec):
        self.seq = genome_rec.seq
        self.desc = genome_rec.description
        self.name = genome_rec.id
        self.length = len(genome_rec.seq)

    @staticmethod
    def parse_genomes(fasta_file):
        genomes_dict = OrderedDict()

        if os.path.exists(fasta_file):
            for genome_rec in SeqIO.parse(fasta_file, "fasta"):
                genomes_dict[genome_rec.id] = Genome(genome_rec)
        else:
            logging.error("Filename {} does not exist.".format(fasta_file))
            raise FileNotFoundError("{} file does not exist.".format(fasta_file))

        return genomes_dict
