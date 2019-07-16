from ..common.RangesHits import RangesHits


def find_overlaps(query, subject, ignore_strand=False, min_overlap=None, min_query_overlap_percentage=None, min_subject_overlap_percentage=None):

    unbinded_query = set(list(range(len(query))))
    unbinded_subject = set(list(range(len(subject))))
    match = []
    query_no_overlap = []
    query_low_overlap = []

    # For each query record
    for qi, q in enumerate(query):
        # For each subject record
        q_no_overlap = True
        q_low_overlap = True
        for si, s in enumerate(subject):
            # Test for overlap
            overlap = q.get_overlap_length(s, ignore_strand)
            filter = False
            if min_overlap is not None:
                filter = filter or overlap < min_overlap
            if min_query_overlap_percentage is not None:
                filter = filter or (overlap/q.width)*100 < min_query_overlap_percentage
            if min_query_overlap_percentage is not None:
                filter = filter or (overlap/s.width)*100 < min_subject_overlap_percentage
            if overlap > 0:
                q_no_overlap = False
            if not filter:
                q_low_overlap = False
                match.append((qi, si))
                unbinded_query.discard(qi)
                unbinded_subject.discard(si)
        query_no_overlap.append(q_no_overlap)
        query_low_overlap.append(q_low_overlap)

    return RangesHits(len(query), len(subject), match, unbinded_query, unbinded_subject, query_no_overlap, query_low_overlap)


def get_unique_ranges(ranges, min_overlap=1):

    hits = find_overlaps(ranges, ranges, ignore_strand=False, min_overlap=min_overlap)

    unique_ranges = set(list(range(len(ranges))))
    for i, j in hits.match:
        if i in unique_ranges:
            if i != j:
                unique_ranges.discard(j)

    return list(hits.unbinded_query) + list(unique_ranges)
