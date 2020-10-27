import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
HMM_MODELS = os.path.join(ROOT_DIR, "data", "models", "hmm", "hmm_all_subfams.hmm")
OUTLIERS = os.path.join(ROOT_DIR, "data", "models", "fasta", "outliers.fasta")
ISFINDER_ORF_DB = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_orf")
ISFINDER_IS_DB = os.path.join(ROOT_DIR, "data", "Blast_db", "isfinder_is")
CONTEXT_SIZE_ORF = 1600
CONTEXT_SIZE_IS = 14000
MIN_GB_OVERLAP = 100
MAX_MERGE_DISTANCE = 700
MIN_HIT_LENGTH = 150
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

INTRAFAMILY_ORF_SIM_THRESHOLD = 0.45
INTERFAMILY_ORF_SIM_THRESHOLD = 0.25

INTRAFAMILY_DNA_SIM_THRESHOLD = 0.70
INTERFAMILY_DNA_SIM_THRESHOLD = 0.50

NUM_THREADS = 4
MIN_EVAL_OVERLAP = 100

FASTA_EXTENSIONS = [".fasta", ".fna", ".ffn"]
GENBANK_EXTENSIONS = [".gb", ".genbank", ".gbff"]

IS_GB_KEYWORDS = ['transposase', 'insertion element', 'mobile element', 'transposon', 'transposable element', 'DDE', 'resolvase', 'recombinase', 'recombination/resolution']

IS_FAMILIES_NAMES = ['IS1', 'IS110', 'IS1182', 'IS1380', 'IS1595', 'IS1634', 'IS21', 'IS256',
                     'IS3', 'IS30', 'IS4', 'IS481', 'IS5', 'IS6', 'IS66', 'IS607', 'IS630', 'IS701', 'IS91',
                     'IS982', 'ISAs1', 'ISAzo13', 'ISH3', 'ISH6', 'ISKra4', 'ISL3', 'ISLre2', 'Tn3', 'ISNCY']

IS_SUBFAMILIES_NAMES = ['ISMhu11', 'IS2', 'IS51', 'IS150', 'IS407', 'IS4', 'IS4Sa', 'IS10', 'IS231',
                        'IS4Sa', 'IS50', 'ISH8', 'ISPepr1', 'IS5', 'IS427', 'IS903', 'IS1031', 'ISH1', 'ISL2',
                        'ISBst12', 'IS1111', 'IS200', 'IS605', 'IS1249', 'ISC1250', 'ISAba11', 'IS942', 'IS1016',
                        'ISH4', 'ISNha5', 'ISNwi1', 'ISPna2', 'ISSod11', 'ISAzba1', 'ISMich2' ]

HYPOTHETICAL_GB_KEYWORDS = ['hypothetical protein', 'predicted protein', 'unknown']

NEUTRAL_GB_KEYWORDS = ['dispersed repetitive unit', 'Tn-like element']
