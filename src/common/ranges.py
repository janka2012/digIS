from ..common.RangesHits import RangesHits


def find_overlaps(query, subject, min_overlap, ignore_strand):

    unbinded_query = set(list(range(len(query))))
    unbinded_subject = set(list(range(len(subject))))
    match = []
    query_max_overlap = []

    # For each query record
    for qi, q in enumerate(query):
        # For each subject record
        max_overlap = 0
        for si, s in enumerate(subject):
            # Test for overlap
            overlap = q.get_overlap_length(s, ignore_strand)
            if overlap > max_overlap:
                max_overlap = overlap
            if overlap >= min_overlap:
                match.append((qi, si))
                unbinded_query.discard(qi)
                unbinded_subject.discard(si)
        query_max_overlap.append(max_overlap)

    return RangesHits(len(query), len(subject), match, unbinded_query, unbinded_subject, query_max_overlap)

