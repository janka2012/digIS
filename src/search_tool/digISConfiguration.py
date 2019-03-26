import os
from definitions import ROOT_DIR


class digISConfiguration:
    def __init__(self, context_size_orf, context_size_is, max_merge_distance, min_gb_overlap):
        self.models = os.path.join(ROOT_DIR, "data", "models", "hmm", "hmm_all_fams.hmm")
        self.isfinder_orf_db = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_orf")
        self.isfinder_is_db = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_is")
        self.context_size_orf = context_size_orf
        self.context_size_is = context_size_is
        self.max_merge_distance = max_merge_distance
        self.min_gb_overlap = min_gb_overlap

