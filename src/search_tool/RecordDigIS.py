from ..common.grange import Grange
from ..search_tool.digISClassification import digISClassification


class RecordDigIS(Grange):
    def __init__(self, genome, chrom, seq_file, genome_len, qid, sid, qstart, qend, start, end, strand, acc):
        super().__init__(genome, chrom, start, end, strand, seq_file, genome_len)
        self.qid = qid
        self.sid = sid
        self.qstart = qstart
        self.qend = qend
        self.acc = acc
        self.classification = digISClassification()

    @classmethod
    def from_csv(cls, csv, genome, chrom, seq_file, genome_len):
        qid = csv['qid']
        sid = csv['sid']
        qstart = int(csv['qstart'])
        qend = int(csv['qend'])
        start = int(csv['sstart'])
        end = int(csv['send'])
        strand = csv['strand']
        acc = csv['acc']
        return cls(genome, chrom, seq_file, genome_len, qid, sid, qstart, qend, start, end, strand, acc)

    @classmethod
    def from_hmmer(cls, hsp, sid, start, end, strand, genome, chrom, seq_file, seq_len):
        return cls(genome, chrom, seq_file, seq_len, hsp.qid, sid, hsp.qstart, hsp.qend, start, end, strand, hsp.acc)

    def should_be_merged(self, other, merge_distance):
        if self.qid == other.qid \
            and not self.has_overlap(other) \
                and self.has_overlap(other, flank=merge_distance):
            return True
        else:
            return False

    def merge(self, other):
        if self.strand != other.strand and self.sid != other.sid:
            raise ValueError('RecordDigIS.merge(): Records can not be merged')

        new_start = min(self.start, other.start)
        new_end = max(self.end, other.end)
        new_len = new_end - new_start + 1

        intersection_length = self.get_overlap_length(other)
        if self.acc > other.acc:
            new_acc = (len(self)*self.acc + (len(other)-intersection_length)*other.acc) / new_len
        else:
            new_acc = ((len(self)-intersection_length)*self.acc + len(other)*other.acc) / new_len

        self.start = new_start
        self.end = new_end
        self.qstart = min(self.qstart, other.qstart)
        self.qend = max(self.qend, other.qend)
        self.qid = '_'.join(list(set(self.qid.split('_') + other.qid.split('_'))))
        self.acc = new_acc

    def __str__(self):
        return "{}, {}, {}, {}, {}, {}, {}, {}, {}".format(self.qid, self.sid, self.qstart, self.qend, self.start,
                                                           self.end, self.strand, self.acc, self.classification)
