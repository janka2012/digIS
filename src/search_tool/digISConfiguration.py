from pathlib import Path
from ..common.misc import init_output_dir


class digISConfiguration:
    def __init__(self, genome_file, context_size_orf, context_size_is, max_merge_distance,
                 genbank_file, min_gb_overlap, out_format, output_dir):
        self.genome_file = genome_file
        self.models = Path("data/models/hmm/hmm_all_subfams.hmm").resolve()
        self.outliers_fasta = Path("data/models/fasta/outliers.fasta").resolve()
        self.isfinder_orf_db = Path("data/Blast_db/isfinder_orf").resolve()
        self.isfinder_is_db = Path("data/Blast_db/isfinder_is").resolve()
        self.context_size_orf = context_size_orf
        self.context_size_is = context_size_is
        self.max_merge_distance = max_merge_distance
        self.genbank_file = genbank_file
        self.min_gb_overlap = min_gb_overlap
        self.out_format = out_format
        self.output_dir = output_dir
        init_output_dir(self.output_dir)

        if self.context_size_orf > self.context_size_is:
            print("Context size ORF is greater than context_size_is. Should be smaller or equal.")
            print("Context size ORF value: {}".format(self.context_size_orf))
            print("Context size IS value: {}".format(self.context_size_is))
            exit(1)

    def __str__(self):
        return "genome file: {}, models: {}, outliers: {}, isfinder orf db: {}, isfinder is db: {}, context orf size: {}, " \
               "context size is: {}, max merge distance: {}, genbank file: {}, min gb overlap {}, outfmt: {}, outdir: {}".format(
              self.genome_file, self.models, self.outliers_fasta, self.isfinder_orf_db, self.isfinder_is_db, self.context_size_orf,
              self.context_size_is, self.max_merge_distance, self.genbank_file, self.min_gb_overlap, self.out_format, self.output_dir)


