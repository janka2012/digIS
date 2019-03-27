from ..blast.BlastHit import BlastHit
from ..blast.BlastHsp import BlastHsp
from ..blast.BlastHspFlat import BlastHspFlat


class BlastQuery:
    def __init__(self, query_id, query_length, hits):
        self.query_id = query_id
        self.query_length = query_length
        self.hits = hits

    @classmethod
    def from_rec(cls, rec):
        hits = []
        for hit in rec.alignments:
            hits.append(BlastHit.from_rec(hit))
        if rec.application == "BLASTX":
            query_length = int(rec.query_length / 3)
        else:
            query_length = rec.query_length
        return cls(rec.query, query_length, hits)

    def get_best_hit(self, query_range=(0, 0), min_overlap=1):
        bhsp = BlastHsp()
        bhspflat = BlastHspFlat()
        for hit in self.hits:
            hsp = hit.get_best_hsp(query_range, min_overlap)
            if bhsp < hsp:
                bhsp = hsp
                bhspflat.set_from_hsp(hsp, self.query_id, self.query_length,
                                      hit.subject_id, hit.subject_length)
        return bhspflat

    def __str__(self):
        hits = '\n'.join(list(str(i) for i in self.hits))
        return "{}, {}\n".format(self.query_id, self.query_length) + hits
