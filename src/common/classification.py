from digIS.src.blast.Blast import Blast
from digIS.src.blast.BlastX import BlastX
from digIS.src.blast.BlastN import BlastN
from digIS.src.common.Classifier import Classifier
from digIS.src.common.ranges import find_overlaps
from digIS.src.common.genbank import read_gb
from digIS.src.genbank.RecordGenbank import RecordGenbank

def get_recs_genbank_overlap_map(self, genbank_records=None, ignore_strand=False):
    hits = find_overlaps(genbank_records, self.recs, self.config.min_gb_overlap, ignore_strand)
    subj_map = hits.get_subject_map()
    recs_genbank = []

    for map in subj_map:
        recs_genbank.append([genbank_records[idx] for idx in map])

    return recs_genbank

def classification(recs, genome, config):
    ds_genbank_recs = [None] * len(recs)
    classification_recs = []

    if recs:
        orf_blast_hits = Blast.get_best_blast_hits_in_range(recs, BlastX,
                                                            config.context_size_orf,
                                                            config.isfinder_orf_db)
        is_dna_blast_hits = Blast.get_best_blast_hits_in_range(recs, BlastN,
                                                               config.context_size_is,
                                                               config.isfinder_is_db)

        if config.genbank_file:
            genbank_recs = list(RecordGenbank(i, genome.name, "chr", genome.file, genome.length)
                                for i in read_gb(config.genbank_file))
            genbank_recs = list(rec for rec in genbank_recs if rec.type not in ['source'])
            ds_genbank_recs = get_recs_genbank_overlap_map(genbank_recs, ignore_strand=True)

        for digis_hit, gb_rec, orf_blast_hit, is_blast_hit in zip(recs, ds_genbank_recs, orf_blast_hits,
                                                                  is_dna_blast_hits):
            classifier = Classifier(digis_hit, gb_rec, orf_blast_hit, is_blast_hit)
            classifier.classify()
            classification_recs.append(classifier)

    return classification_recs

