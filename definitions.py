import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
HMM_MODELS = os.path.join(ROOT_DIR, "data", "models", "hmm", "hmm_all_subfams.hmm")
OUTLIERS = os.path.join(ROOT_DIR, "data", "models", "fasta", "outliers.fasta")
ISFINDER_ORF_DB = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_orf_all")
ISFINDER_IS_DB = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_is")
CONTEXT_SIZE_ORF = 1600
CONTEXT_SIZE_IS = 14000
MIN_GB_OVERLAP = 100
MAX_MERGE_DISTANCE = 700  # nt
MIN_HIT_LENGTH = 150  # nt
CURRATED_CUTOFF = False
OUTLIERS_EVALUE = 0.001

BLASTN_GAPOPEN = 5
BLASTN_GAPEXTEND = 2
BLASTN_WORDSIZE = 11
BLASTN_EVALUE = 0.001

BLASTX_GAPOPEN = 11
BLASTX_GAPEXTEND = 1
BLASTX_WORDSIZE = 3
BLASTX_EVALUE = 0.001

NUM_THREADS = 0

IS_GB_KEYWORDS = ['transposase', 'resolvase', 'recombinase', 'recombination/resolution',
                  'insertion element', 'mobile element', 'transposon', 'transposable element', 'DDE']

IS_FAMILIES_NAMES = ['IS1', 'IS110', 'IS1182', 'IS1380', 'IS1595', 'IS1634', 'IS200', 'IS605', 'IS21', 'IS256',
                     'IS3', 'IS30', 'IS4', 'IS481', 'IS5', 'IS6', 'IS66', 'IS607', 'IS630', 'IS701', 'IS91',
                     'IS982', 'ISAs1', 'ISAzo13', 'ISH3', 'ISH6', 'ISKra4', 'ISL3', 'ISLe2', 'Tn3', 'ISNCY']

HYPOTHETICAL_GB_KEYWORDS = ['hypothetical protein', 'predicted protein', 'unknown']