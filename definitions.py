import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
HMM_MODELS = os.path.join(ROOT_DIR, "data", "models", "hmm", "hmm_all_subfams.hmm")
OUTLIERS = os.path.join(ROOT_DIR, "data", "models", "fasta", "outliers.fasta")
ISFINDER_ORF_DB = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_orf")
ISFINDER_IS_DB = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_is")

