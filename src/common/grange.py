from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging
import copy


class Grange:
    def __init__(self, genome, chrom, start, end, strand, seq_file, genome_len):
        self.genome = genome
        self.chr = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.seq_file = seq_file
        self.genome_len = genome_len
        self.width = self.__len__()

    def get_flank_range(self, flank):
        if flank + flank + self.width > self.genome_len:
            logging.error("Element with flank regions is wider than genome length ({}, {}, {}, {}, {}). ".format(self.start, self.end, self.genome_len, flank, self.width) + self.__str__())
            exit(1)
        flank_start = self.start - flank
        if flank_start <= 0:
            flank_start += self.genome_len
        flank_end = self.end + flank
        if flank_end > self.genome_len:
            flank_end -= self.genome_len
        return Grange(self.genome, self.chr, flank_start, flank_end, self.strand, self.seq_file, self.genome_len)

    def shift_left(self, size):
        self.start -= size
        if self.start <= 0:
            self.start += self.genome_len
        self.end -= size
        if self.end <= 0:
            self.end += self.genome_len

    def has_overlap(self, other, ignore_strand=False, flank=0):
        return self.get_overlap_length(other, ignore_strand, flank) > 0

    def get_overlap_length(self, other, ignore_strand=False, flank=0):
        other_range = copy.copy(other)
        new_range = self.get_flank_range(flank)
        size = min(new_range.start, other_range.start)-1
        new_range.shift_left(size)
        other_range.shift_left(size)
        if new_range.start <= other_range.start:
            overlap = min(new_range.end-other_range.start+1, other.width)
        else:
            overlap = min(other_range.end-new_range.start+1, new_range.width)

        if self.strand != other.strand and not ignore_strand:
            overlap = 0

        # print(self.start, self.end, self.strand, other.start, other.end, other.strand, flank, overlap)
        return overlap

    def is_inside(self, other, ignore_strand=False):
        return self.get_overlap_length(other, ignore_strand, 0) == self.width

    def get_flank_lengths(self, flank):
        return flank, flank

    def get_sequence(self, flank=0, protein=False):
        record = SeqIO.read(self.seq_file, "fasta")
        new_range = self.get_flank_range(flank)
        if new_range.start <= new_range.end:
            seq = record.seq[new_range.start-1:new_range.end]
        else: # element crossing the genome boundary
            seq = record.seq[new_range.start-1:self.genome_len] + record.seq[0:new_range.end]

        if self.strand == '-':
            seq = seq.reverse_complement()
        if protein:
            seq = seq.translate(table=11)

        return SeqRecord(seq, id=record.id, description='')

    def __str__(self):
        return "{} {} {} {} {}".format(self.genome, self.chr, self.start, self.end, self.strand)

    def __len__(self):
        out_len = self.end - self.start + 1
        if out_len <= 0: # element crossing the genome boundary
            out_len += self.genome_len
        return out_len