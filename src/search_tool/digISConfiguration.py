import definitions

from ..common.misc import init_output_dir


class digISConfiguration:
    def __init__(self, genome_file, genbank_file, min_gb_overlap, out_format, output_dir):
        self.genome_file = genome_file
        self.models = definitions.HMM_MODELS
        self.outliers_fasta = definitions.OUTLIERS
        self.isfinder_orf_db = definitions.ISFINDER_ORF_DB
        self.isfinder_is_db = definitions.ISFINDER_IS_DB
        self.context_size_orf = 1600
        self.context_size_is = 14000
        self.max_merge_distance = 700
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


