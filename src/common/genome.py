import os

import src.common.sequence as seq

class Genome:

    def __init__(self, genome_rec, output_dir):
        self.seq = genome_rec.seq
        self.desc = genome_rec.description
        self.name = genome_rec.id
        self.length = len(genome_rec.seq)
        self.orf_db = os.path.join(output_dir, "pep", self.name + ".pep")
        seq.translate_dna_seq_biopython(seqrec=genome_rec, outseq=self.orf_db)
