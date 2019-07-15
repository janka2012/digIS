import logging
import os

from collections import Counter, OrderedDict
from Bio import SeqIO
from ..common.sequence import translate_dna_seq_biopython


class Genome:

    def __init__(self, genome_rec, translate=True, output_dir=None):
        self.seq = genome_rec.seq
        self.desc = genome_rec.description
        self.name = genome_rec.id
        self.length = len(genome_rec.seq)
        self.orf_db = os.path.join(output_dir, "pep", self.name + ".pep")
        if translate:
            translate_dna_seq_biopython(seqrec=genome_rec, outseq=self.orf_db)

    @staticmethod
    def parse_genomes(fasta_file, translate=True, output_dir=None):
        genomes_dict = OrderedDict()
        genome_ids = []
        if os.path.exists(fasta_file):
            for genome_rec in SeqIO.parse(fasta_file, "fasta"):
                genome_ids.append(genome_rec.id)
                genomes_dict[genome_rec.id] = Genome(genome_rec, translate, output_dir)

            # check if there are duplicates in genome_ids
            duplicated_ids = [item for item, count in Counter(genome_ids).items() if count > 1]

            if duplicated_ids:
                raise ValueError("Multifasta file contains following duplicated ids: {}.\n"
                                 "Ids have to be unique.".format(duplicated_ids))
            print(duplicated_ids)

        else:
            logging.error("Filename {} does not exist.".format(fasta_file))
            raise FileNotFoundError("{} file does not exist.".format(fasta_file))

        return genomes_dict
