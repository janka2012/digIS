from digIS.src.blast.Blast import Blast
from digIS.src.blast.BlastX import BlastX
from digIS.src.blast.BlastN import BlastN
from digIS.src.common.Classifier import Classifier
from digIS.src.common.ranges import find_overlaps

def get_recs_genbank_overlap_map(recs, genbank_records=None, min_overlap=100, ignore_strand=False):
    hits = find_overlaps(genbank_records, recs, min_overlap, ignore_strand)
    subj_map = hits.get_subject_map()
    recs_genbank = []

    for map_list in subj_map:
        recs_genbank.append([genbank_records[idx] for idx in map_list])

    return recs_genbank

def classification(recs, gb_recs, context_size_orf, context_size_is, min_overlap, isfinder_orf_db, isfinder_is_db):
    ds_genbank_recs = [None] * len(recs)
    classification_recs = []

    if recs:
        orf_blast_hits = Blast.get_best_blast_hits_in_range(recs, BlastX, context_size_orf, isfinder_orf_db)
        is_dna_blast_hits = Blast.get_best_blast_hits_in_range(recs, BlastN, context_size_is, isfinder_is_db)
        if gb_recs is not None:
            ds_genbank_recs = get_recs_genbank_overlap_map(recs, gb_recs, min_overlap, ignore_strand=True)

        for digis_hit, gb_rec, orf_blast_hit, is_blast_hit in zip(recs, ds_genbank_recs, orf_blast_hits,
                                                                  is_dna_blast_hits):
            classifier = Classifier(digis_hit, gb_rec, orf_blast_hit, is_blast_hit)
            classifier.classify()
            classification_recs.append(classifier)

    return classification_recs

