from ..blast.BlastHit import BlastHit
from ..blast.BlastHsp import BlastHsp
from ..blast.BlastHspFlat import BlastHspFlat


class BlastQuery:
    def __init__(self, query_id, query_length, hits, application=''):
        self.query_id = query_id
        self.query_length = query_length
        self.hits = hits
        self.application = application

    @classmethod
    def from_rec(cls, rec):
        hits = []
        for hit in rec.alignments:
            hits.append(BlastHit.from_rec(hit))
        query_length = rec.query_length
        return cls(rec.query, query_length, hits, rec.application)

    def get_best_hit(self, query_range=(0, 0), min_overlap=1, positive_subject_strand_only=False):
        bhsp = BlastHsp()
        bhspflat = BlastHspFlat()
        for hit in self.hits:
            hsp = hit.get_best_hsp(query_range, min_overlap, positive_subject_strand_only)
            if bhsp < hsp:
                bhsp = hsp
                bhspflat.set_from_hsp(hsp, self.query_id, self.query_length,
                                      hit.subject_id, hit.subject_length, self.application)
        return bhspflat

    def get_max_hit(self, query_range=(0, 0), min_overlap=1, positive_subject_strand_only=False):
        query_start, query_end = 0,0

        for hit in self.hits:
            hit_query_start, hit_query_end = hit.get_max_hsp(query_range, min_overlap, positive_subject_strand_only)
            if hit_query_start == 0 and hit_query_end == 0:
                continue
            if query_start == 0 and query_end == 0:
                query_start, query_end = hit_query_start, hit_query_end
            else:
                query_start = min(query_start, hit_query_start)
                query_end = max(query_end, hit_query_end)

        return query_start, query_end


    def __str__(self):
        hits = '\n'.join(list(str(i) for i in self.hits))
        return "{}, {}\n".format(self.query_id, self.query_length) + hits
