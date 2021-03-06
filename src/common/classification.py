from ..blast.Blast import Blast
from ..blast.BlastX import BlastX
from ..blast.BlastN import BlastN
from ..common.Classifier import Classifier
from ..common.ranges import find_overlaps


def get_recs_genbank_overlap_map(recs, genbank_records=None, min_gb_overlap=None, ignore_strand=False):
    hits = find_overlaps(genbank_records, recs, ignore_strand=ignore_strand, min_overlap=min_gb_overlap)
    subj_map = hits.get_subject_map()
    recs_genbank = []

    for map_list in subj_map:
        recs_genbank.append([genbank_records[idx] for idx in map_list])

    return recs_genbank


def classification(recs, gb_recs, context_size_orf, context_size_is, min_gb_overlap, isfinder_orf_db, isfinder_is_db):
    ds_genbank_recs = [None] * len(recs)
    classification_recs = []

    if recs:
        orf_blast_hits = Blast.get_best_blast_hits_in_range(recs, BlastX, context_size_orf, isfinder_orf_db, min_overlap=min_gb_overlap, positive_subject_strand_only=True)
        is_dna_blast_hits = Blast.get_best_blast_hits_in_range(recs, BlastN, context_size_is, isfinder_is_db, min_overlap=min_gb_overlap, positive_subject_strand_only=True)
        if gb_recs is not None:
            ds_genbank_recs = get_recs_genbank_overlap_map(recs, gb_recs, min_gb_overlap, ignore_strand=True)

        for digis_hit, gb_rec, orf_blast_hit, is_blast_hit in zip(recs, ds_genbank_recs, orf_blast_hits,
                                                                  is_dna_blast_hits):
            classifier = Classifier(digis_hit, gb_rec, orf_blast_hit, is_blast_hit)
            classifier.classify()
            classification_recs.append(classifier)

    return classification_recs

