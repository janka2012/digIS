import os
from ..common.sequence import get_seqlen
from ..common.sequence import translate_dna_seq_biopython


class Genome:

    def __init__(self, genome_file, output_dir):
        self.file = genome_file
        self.name = os.path.splitext(os.path.basename(self.file))[0]
        self.length = get_seqlen(self.file)
        self.orf_db = os.path.join(output_dir, "pep", self.name + ".pep")
        translate_dna_seq_biopython(sequence=self.file, outseq=self.orf_db)
