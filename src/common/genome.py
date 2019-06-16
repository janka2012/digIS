import os

from Bio import SeqIO
from ..common.sequence import translate_dna_seq_biopython


class Genome:

    def __init__(self, genome_rec, output_dir):
        self.seq = genome_rec.seq
        self.desc = genome_rec.description
        self.name = genome_rec.id
        self.length = len(genome_rec.seq)
        self.orf_db = os.path.join(output_dir, "pep", self.name + ".pep")
        translate_dna_seq_biopython(seqrec=genome_rec, outseq=self.orf_db)

    @staticmethod
    def parse_genomes(fasta_file, output_dir):
        genomes_dict = {}
        if os.path.exists(fasta_file):
            for genome in SeqIO.parse(fasta_file, "fasta"):
                # TODO merge output genome name with filename?
                genomes_dict[genome.id] = Genome(genome, output_dir)

        return genomes_dict

