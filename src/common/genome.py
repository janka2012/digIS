import os

from Bio import SeqIO
from ..common.sequence import translate_dna_seq_biopython


class Genome:

    def __init__(self, genome_rec, translate=True, output_dir=None):
        self.seq = genome_rec.seq
        self.desc = genome_rec.description
        self.name = genome_rec.id
        self.length = len(genome_rec.seq)
        if translate:
            self.orf_db = os.path.join(output_dir, "pep", self.name + ".pep")
            translate_dna_seq_biopython(seqrec=genome_rec, outseq=self.orf_db)

    @staticmethod
    def parse_genomes(fasta_file, translate=True, output_dir=None):
        genomes_dict = {}
        if os.path.exists(fasta_file):
            for genome in SeqIO.parse(fasta_file, "fasta"):
                genomes_dict[genome.id] = Genome(genome, translate, output_dir)

        return genomes_dict

