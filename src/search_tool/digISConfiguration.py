import os

from ..common.misc import init_output_dir
from definitions import ROOT_DIR


class digISConfiguration:
    def __init__(self, genome_file, context_size_orf, context_size_is, max_merge_distance,
                 genbank_file, min_gb_overlap, out_format, output_dir):
        self.genome_file = genome_file
        self.models = os.path.join(ROOT_DIR, "data", "models", "hmm", "hmm_all_fams.hmm")
        self.outliers_fasta = os.path.join(ROOT_DIR, "data", "models", "fasta", "outliers.fasta")
        self.isfinder_orf_db = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_orf")
        self.isfinder_is_db = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_is")
        self.context_size_orf = context_size_orf
        self.context_size_is = context_size_is
        self.max_merge_distance = max_merge_distance
        self.genbank_file = genbank_file
        self.min_gb_overlap = min_gb_overlap
        self.out_format = out_format
        self.output_dir = output_dir
        init_output_dir(self.output_dir)

