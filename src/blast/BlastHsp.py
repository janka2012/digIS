class BlastHsp:
    def __init__(self, query_start=0, query_end=0, subject_start=0, subject_end=0,
                 subject_strand='+', score=0.0, identities=0, positives=0):
        self.query_start = query_start
        self.query_end = query_end
        self.subject_start = subject_start
        self.subject_end = subject_end
        self.subject_strand = subject_strand
        self.score = score
        self.identities = identities
        self.positives = positives

    @classmethod
    def from_rec(cls, hsp):
        if hsp.sbjct_start < hsp.sbjct_end:
            subject_strand = '+'
            subject_start, subject_end = hsp.sbjct_start, hsp.sbjct_end
        else:
            subject_strand = '-'
            subject_end, subject_start = hsp.sbjct_start, hsp.sbjct_end
        return cls(hsp.query_start, hsp.query_end, subject_start, subject_end,
                   subject_strand, hsp.score, hsp.identities, hsp.positives)

    def __lt__(self, other):
        return self.score < other.score

    def __eq__(self, other):
        return self.score == other.score

    def __add__(self, other):
        return BlastHsp(self.query_start, other.query_end, self.subject_start, other.subject_end, self.subject_strand,
                        self.score + other.score, self.identities + other.identities, self.positives + other.positives)

    def __str__(self):
        return "Sc: {:1.1f}, Id: {}, Pos: {}, QStart: {}, QEnd: {}, SStart: {}, SEnd: {}, SStrand: {}"\
            .format(self.score, self.identities, self.positives, self.query_start,
                    self.query_end, self.subject_start, self.subject_end, self.subject_strand)
