from copy import copy

from ..blast.BlastHsp import BlastHsp


class BlastHit:
    def __init__(self, subject_id, subject_length, hsps):
        self.subject_id = subject_id
        self.subject_length = subject_length
        self.hsps = hsps

    @classmethod
    def from_rec(cls, rec):
        subject_id = rec.hit_id
        subject_length = rec.length
        hsps = []
        for hsp in rec.hsps:
            hsps.append(BlastHsp.from_rec(hsp))
        return cls(subject_id, subject_length, hsps)

    def merge_hsps(self, max_space):
        hsps = copy(self.hsps)
        change = True
        merged_idx = []
        while change:
            change = False
            merged = []
            for i, hi in enumerate(hsps):
                for j, hj in enumerate(hsps):
                    if i != j and not (i, j) in merged_idx \
                            and hi.query_end <= hj.query_start \
                            and hi.subject_end <= hj.subject_start \
                            and hi.subject_strand == hj.subject_strand \
                            and hj.query_start - hi.query_end + 1 <= max_space \
                            and hj.subject_start - hi.subject_end + 1 <= max_space:
                        merged.append(hi+hj)
                        change = True
                        merged_idx.append((i, j))
            hsps = hsps + merged
        return hsps

    def get_best_hsp(self, query_range=(0, 0), min_overlap=1, positive_subject_strand_only=False):
        bhsp = BlastHsp()

        for hsp in self.hsps:
            if positive_subject_strand_only and hsp.subject_strand == "-":
                continue
            if query_range == (0, 0):
                if bhsp < hsp:
                    bhsp = hsp
            else:
                max_start = max(hsp.query_start, query_range[0])
                min_end = min(hsp.query_end, query_range[1])
                overlap_len = min_end - max_start + 1

                if bhsp < hsp and overlap_len >= min_overlap:
                    bhsp = hsp
        return bhsp

    def get_max_hsp(self, query_range=(0, 0), min_overlap=1, positive_subject_strand_only=False):
        query_start, query_end = 0,0

        for hsp in self.hsps:
            if positive_subject_strand_only and hsp.subject_strand == "-":
                continue
            max_start = max(hsp.query_start, query_range[0])
            min_end = min(hsp.query_end, query_range[1])
            overlap_len = min_end - max_start + 1

            if query_range == (0, 0) or overlap_len >= min_overlap:
                if query_start == 0 and query_end == 0:
                    query_start = hsp.query_start
                    query_end = hsp.query_end
                else:
                    query_start = min(query_start, hsp.query_start)
                    query_end = max(query_end, hsp.query_end)

        return query_start, query_end

    def __str__(self):
        hsps = '\n'.join(list(str(i) for i in self.hsps))
        return "{}, {}\n".format(self.subject_id, self.subject_length) + hsps
