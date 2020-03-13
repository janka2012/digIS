from ..common.grange import Grange


class RecordDigIS(Grange):
    def __init__(self, genome_name, chrom, genome_seq, genome_len, qid, sid, qstart, qend, start, end, strand, acc, score, evalue):
        self.qid = qid
        self.sid = sid
        self.qstart = qstart
        self.qend = qend
        self.acc = acc
        self.score = score
        self.evalue = evalue
        super().__init__(genome_name, chrom, start, end, strand, genome_seq, genome_len)

    @classmethod
    def from_csv(cls, csv, genome_name, chrom, genome_seq, genome_len):
        qid = csv['qid']
        sid = csv['sid']
        qstart = int(csv['qstart'])
        qend = int(csv['qend'])
        start = int(csv['sstart'])
        end = int(csv['send'])
        strand = csv['strand']
        acc = csv['acc']
        score = csv['score']
        evalue = csv['evalue']
        return cls(genome_name, chrom, genome_seq, genome_len, qid, sid, qstart, qend, start, end, strand, acc, score, evalue)

    @classmethod
    def from_hmmer(cls, hsp, sid, start, end, strand, genome_name, chrom, genome_seq, seq_len):
        return cls(genome_name, chrom, genome_seq, seq_len, hsp.qid, sid, hsp.qstart, hsp.qend, start, end, strand, float(hsp.acc), float(hsp.dom_bitscore), float(hsp.dom_evalue))

    # Regurements for merge in distance
    # - the same strand
    # - the same query_id (hmm model/outlier)
    # - continuous fragments with respect to model
    def should_be_merged_distance(self, other, merge_distance):
        
        continuous_fragments = (self.strand == '+' and self.start < other.start and self.qend <= other.qstart) or \
            (self.strand == '+' and other.start < self.start and other.qend <= self.qstart) or \
                (self.strand == '-' and other.start < self.start and self.qend <= other.qstart) or \
                    (self.strand == '-' and self.start < other.start and other.qend <= self.qstart)
        
        if self.qid == other.qid and not self.has_overlap(other) \
                and self.has_overlap(other, flank=merge_distance) and continuous_fragments:
            return True
        else:
            return False

    def merge(self, other, merge_type):
        if self.strand != other.strand or self.sid != other.sid:
            raise ValueError('RecordDigIS.merge(): Records can not be merged')

        # if both hits from the same

        new_start = min(self.start, other.start)
        new_end = max(self.end, other.end)
        new_len = new_end - new_start + 1

        intersection_length = self.get_overlap_length(other)
        if self.acc > other.acc:
            new_acc = (len(self)*self.acc + (len(other)-intersection_length)*other.acc) / new_len
        else:
            new_acc = ((len(self)-intersection_length)*self.acc + len(other)*other.acc) / new_len

        if merge_type == "distance":
            new_score = self.score + other.score
        elif merge_type == "overlap":
            new_score = max(self.score, other.score)

        self.start = new_start
        self.end = new_end
        self.qstart = min(self.qstart, other.qstart)
        self.qend = max(self.qend, other.qend)
        self.qid = '-'.join(list(set(self.qid.split('-') + other.qid.split('-'))))
        self.acc = new_acc
        self.score = new_score
        self.evalue = min(self.evalue, other.evalue)

    @classmethod
    def get_csv_header(cls):
        return ["qid", "qstart", "qend", "sid", "sstart", "send", "strand", "acc", "score", "evalue"]

    def to_csv(self):
        return [self.qid, self.qstart, self.qend, self.sid, self.start, self.end, self.strand, round(self.acc, 2), round(self.score, 2), self.evalue]

    def __str__(self):
        return "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(self.qid, self.qstart, self.qend,
                                                       self.sid, self.start,  self.end,
                                                       self.strand, self.acc, self.score, self.evalue)
