class HmmerHsp:

    def __init__(self, hsp, query_len):
        self.qid = hsp.query_id
        self.sid = hsp.hit_id
        self.qstart = hsp.query_start + 1
        self.qend = hsp.query_end
        self.sstart = hsp.hit_start + 1
        self.send = hsp.hit_end
        self.sstart_env = hsp.env_start + 1
        self.send_env = hsp.env_end
        self.acc = float(hsp.acc_avg)
        self.evalue = float(hsp.evalue)
        self.bitscore = float(hsp.bitscore)
        self.query_coverage = (self.qend - self.qstart + 1) / query_len

    def __str__(self):
        return "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(self.qid, self.sid, self.qstart, self.qend, self.sstart,
                                                               self.send, self.acc, self.evalue, self.bitscore,
                                                               self.query_coverage)

    def __lt__(self, other):
        return self.bitscore < other.bitscore
