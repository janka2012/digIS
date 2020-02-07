import definitions

from ..common.misc import init_output_dir

class digISConfiguration:
    def __init__(self, genome_file, genbank_file, output_dir, currated_cutoff=None, outliers_evalue=None):
        self.genome_file = genome_file
        self.models = definitions.HMM_MODELS
        self.outliers_fasta = definitions.OUTLIERS
        self.isfinder_orf_db = definitions.ISFINDER_ORF_DB
        self.isfinder_is_db = definitions.ISFINDER_IS_DB
        self.context_size_orf = definitions.CONTEXT_SIZE_ORF
        self.context_size_is = definitions.CONTEXT_SIZE_IS
        self.max_merge_distance = definitions.MAX_MERGE_DISTANCE
        self.min_hit_length = definitions.MIN_HIT_LENGTH
        self.min_gb_overlap = definitions.MIN_GB_OVERLAP
        self.currated_cutoff = currated_cutoff if currated_cutoff else definitions.CURRATED_CUTOFF
        self.outliers_evalue = outliers_evalue if outliers_evalue else definitions.OUTLIERS_EVALUE
        self.genbank_file = genbank_file
        self.output_dir = output_dir
        init_output_dir(self.output_dir)

        if self.context_size_orf > self.context_size_is:
            msg = "Context size ORF is greater than context_size_is. Should be smaller or equal.\n"
            msg += "Context size ORF value: {}\n".format(self.context_size_orf)
            msg += "Context size IS value: {}\n".format(self.context_size_is)

            raise ValueError(msg)

    def __str__(self):
        return "genome file: {}, models: {}, outliers: {}, isfinder orf db: {}, isfinder is db: {}, context orf size: {}, " \
               "context size is: {}, max merge distance: {}, genbank file: {}, min gb overlap {}, outdir: {}".format(
              self.genome_file, self.models, self.outliers_fasta, self.isfinder_orf_db, self.isfinder_is_db, self.context_size_orf,
              self.context_size_is, self.max_merge_distance, self.genbank_file, self.min_gb_overlap, self.output_dir)


