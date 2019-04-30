class BlastHspFlat:

    def __init__(self):
        self.score = 0.0
        self.identities = 0
        self.positives = 0
        self.query_id = ""
        self.subject_id = ""
        self.query_start = 0
        self.query_end = 0
        self.subject_start = 0
        self.subject_end = 0
        self.subject_strand = '+'
        self.query_len = 0
        self.subject_len = 0
        self.query_identity = 0.0
        self.subject_identity = 0.0
        self.query_coverage = 0.0
        self.subject_coverage = 0.0
        self.shorter_identity = 0.0
        self.application = ''

    def set_from_hsp(self, hsp, query_id, query_len, subject_id, subject_len, application):
        self.score = hsp.score
        self.identities = hsp.identities
        self.positives = hsp.positives
        self.query_start = hsp.query_start
        self.query_end = hsp.query_end
        self.subject_start = hsp.subject_start
        self.subject_end = hsp.subject_end
        self.subject_strand = hsp.subject_strand
        self.query_id = query_id
        self.query_len = query_len
        self.subject_id = subject_id
        self.subject_len = subject_len
        self.application = application

        self.query_coverage = (self.query_end - self.query_start + 1) / self.query_len
        self.subject_coverage = (self.subject_end - self.subject_start + 1) / self.subject_len
        self.subject_identity = self.identities / self.subject_len

        if self.application == "BLASTX":
            self.query_identity = self.identities / int(self.query_len/3)
            self.shorter_identity = self.identities / min(int(self.query_len/3), self.subject_len)
        else:
            self.query_identity = self.identities / self.query_len
            self.shorter_identity = self.identities / min(self.query_len, self.subject_len)

    def __lt__(self, other):
        return self.score < other.score

    def __eq__(self, other):
        return self.score == other.score

    def __str__(self):
        return 'Qlen: {:4}, DBlen: {:4}, Score: {:6.1f}, Pos: {:4}, Id: {:4}, Qid: {:1.2f}, ' \
               'DBid: {:1.2f}, Qcov: {:1.2f}, DBcov: {:1.2f}, Qid: {}, DBid: {}'.format(
            self.query_len, self.subject_len, self.score, self.positives, self.identities, self.query_identity,
            self.subject_identity, self.query_coverage, self.subject_coverage, self.query_id, self.subject_id)
