class RangesHits:

    def __init__(self, queries, subjects, match, unbinded_query, unbinded_subject, query_no_overlap, query_low_overlap):
        self.queries = queries
        self.subjects = subjects
        self.match = match
        self.unbinded_query = unbinded_query
        self.unbinded_subject = unbinded_subject
        self.query_no_overlap = query_no_overlap
        self.query_low_overlap = query_low_overlap

    def get_subject_map(self):
        subject_map = [[]] * self.subjects
        subject_match = set(subject for query, subject in self.match)
        for k in subject_match:
            for i, j in self.match:
                if j == k:
                    subject_map[j] = subject_map[j] + [i]
        return subject_map
