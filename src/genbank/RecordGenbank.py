from src.common.grange import Grange


class RecordGenbank(Grange):

    def __init__(self, rec, genome="", chrom="chr", seq_file="", genome_len=0):
        strand = "+" if rec.location.strand == 1 else "-"
        if rec.location.strand == 1:
            strand = "+"
            start_pos = rec.location.parts[0].start
            end_pos = rec.location.parts[len(rec.location.parts) - 1].end
        else:
            strand = "-"
            start_pos = rec.location.parts[len(rec.location.parts) - 1].start
            end_pos = rec.location.parts[0].end
        super().__init__(genome, chrom, start_pos+1, end_pos, strand, seq_file, genome_len)
        self.type = rec.type
        self.product = ", ".join(rec.qualifiers['product']) if self.type == "CDS" else ""
        self.qualifiers = rec.qualifiers

    def __str__(self):
        return "{} {} {} {} {} {}".format(self.genome, self.type, self.start, self.end, self.strand, self.qualifiers)
