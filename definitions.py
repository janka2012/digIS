import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
HMM_MODELS = os.path.join(ROOT_DIR, "data", "models", "hmm", "hmm_all_subfams.hmm")
OUTLIERS = os.path.join(ROOT_DIR, "data", "models", "fasta", "outliers.fasta")
ISFINDER_ORF_DB = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_orf")
ISFINDER_IS_DB = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_is")
CONTEXT_SIZE_ORF = 1600
CONTEXT_SIZE_IS = 14000
MIN_GB_OVERLAP = 100
MAX_MERGE_DISTANCE = 700  # nt
MIN_HIT_LENGTH = 150  # nt

