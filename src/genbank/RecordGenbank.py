from ..common.grange import Grange


class RecordGenbank(Grange):

    def __init__(self, rec, genome_name="", chrom="chr", genome_seq="", genome_len=0):
        if rec.location.strand == 1:
            strand = "+"
            start_pos = rec.location.parts[0].start
            end_pos = rec.location.parts[len(rec.location.parts) - 1].end
        else:
            strand = "-"
            start_pos = rec.location.parts[len(rec.location.parts) - 1].start
            end_pos = rec.location.parts[0].end
        super().__init__(genome_name, chrom, start_pos+1, end_pos, strand, genome_seq, genome_len)
        self.type = rec.type
        if 'product' in rec.qualifiers:
            self.product = ", ".join(rec.qualifiers['product']) if self.type == "CDS" else ""
        elif 'note' in rec.qualifiers:
            self.product = ", ".join(rec.qualifiers['note']) if self.type == "CDS" else ""
        else:
            self.product = ""
        self.qualifiers = rec.qualifiers

    def __str__(self):
        return "{} {} {} {} {} {}".format(self.genome_name, self.type, self.start, self.end, self.strand, self.qualifiers)
