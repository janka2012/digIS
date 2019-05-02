import os
import re
import subprocess
import sys

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def transform_range(start, end, frame, seqlen):
    offset = [0, 1, 2, 2, 1, 0][frame-1]
    if frame <= 3:
        start_pos = (start-1)*3 + offset + 1
        end_pos = end*3 + offset
    else:
        start_pos = seqlen - (end*3 + offset) + 1
        end_pos = seqlen - ((start-1)*3 + offset)

    return start_pos, end_pos


def translate_dna_seq_biopython(sequence, outseq):

    with open(outseq, 'w') as aa_fa:
        for dna_record in SeqIO.parse(sequence, 'fasta'):

            aa_seqs = []
            dna_seqs = [dna_record.seq, dna_record.seq.reverse_complement()]
            for dna_seq in dna_seqs:
                for frame in range(1, 4):

                    # correction of end position such that the length is multiple of 3
                    start = frame - 1
                    excess = (len(dna_seq) - start) % 3
                    end = len(dna_seq) - excess
                    seq = dna_seq[start:end]
                    aa_seqs.append(seq.translate(table=11))

            for frame, aa_seq in enumerate(aa_seqs, start=1):
                seq_id = dna_record.id + "_" + str(frame)
                description = seq_id + dna_record.description.replace(dna_record.id, "")
                aa_record = SeqRecord(aa_seq, id=seq_id, name=seq_id, description=description)
                SeqIO.write(aa_record, aa_fa, 'fasta')


def translate_dna_seq(sequence, outseq):

    if not __is_non_zero_file(outseq) or not os.path.isfile(outseq):

        cmd = ["transeq", "-sequence", sequence, "-outseq", outseq, "-frame", "6", "-table", "11"]
        try:
            subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)
        except subprocess.CalledProcessError as e:
            print("An error occurred when calling {}.".format(cmd))
            print(e)


def __is_non_zero_file(filepath):
    return os.path.isfile(filepath) and os.path.getsize(filepath) > 0


def get_seqlen(filename):
    rec = SeqIO.read(filename, "fasta")
    return len(rec.seq)


def get_ids_from_fasta(filename):

    ids = set()
    records = SeqIO.parse(filename, "fasta")
    for rec in records:
        rec_id = rec.id.split("|")[0]
        rec_id = rec_id.split(".")[0]
        ids.add(rec_id)

    return ids


def get_full_ids_from_fasta(filename):

    ids = []
    records = SeqIO.parse(filename, "fasta")
    for rec in records:
        # rec_id = rec.id.split("|")[0]
        # rec_id = rec_id.split(".")[0]
        ids.append(rec.id)

    return ids

def filter_fasta_re(in_file, out_file, regexp):
    out_recs = []
    records = SeqIO.parse(in_file, "fasta")
    for rec in records:
        if re.search(regexp, rec.id, re.IGNORECASE):
            out_recs.append(rec)
    SeqIO.write(out_recs, out_file, "fasta")


def get_sequence_record(filename, start, end, strand, protein=True):
    record = SeqIO.read(filename, "fasta")
    seq = record.seq[start-1:end]

    if strand == '-':
        seq = seq.reverse_complement()
    if protein:
        seq = seq.translate(table=11)

    return SeqRecord(seq, id=record.id, description='')


def get_sequence_record_ids(filename, ids):
    out = []
    ids_set = set(ids)
    recs = SeqIO.parse(filename, "fasta")
    for rec in recs:
        if rec.id in ids_set:
            out.append(rec)

    return out


def get_sequence_record_id(filename, rec_id):
    out = []
    recs = SeqIO.parse(filename, "fasta")
    for rec in recs:
        if rec_id in rec.id:
            out.append(rec)

    return out


def get_sixframe_record(filename, start, end):

    # correction of end position such that the length is multiple of 3
    record = SeqIO.read(filename, "fasta")
    s = start
    e = end + 3 - ((end - start + 1) % 3)

    all_recs = []
    for i in [0, 1, 2]:
        seq = record.seq[s - 1 + i:e + i]
        all_recs.append(SeqRecord(seq.translate(table=11), id=record.id, description='Frame: ' + str(i)))
        rc_seq = seq.reverse_complement()
        all_recs.append(SeqRecord(rc_seq.translate(table=11), id=record.id, description='Frame: ' + str(i + 3)))

    return all_recs


def merge_all_fasta_files(dir_path, out_file):
    all_recs = []
    for (dirpath, dirnames, filenames) in os.walk(dir_path):
        for filename in filenames:
            records = SeqIO.parse(os.path.join(dirpath, filename), "fasta")
            for rec in records:
                all_recs.append(rec)
    SeqIO.write(all_recs, out_file, "fasta")


def merge_fasta_files(filenames, out_file):
    all_recs = []
    for filename in filenames:
        records = SeqIO.parse(filename, "fasta")
        for rec in records:
            all_recs.append(rec)
    SeqIO.write(all_recs, out_file, "fasta")


def get_max_seq_len(filename):
    max_len = 0
    records = SeqIO.parse(filename, "fasta")
    for rec in records:
        max_len = max(max_len, len(rec.seq))

    return max_len


def get_maxlen_seq(filename):
    max_len = 0
    out = None
    records = SeqIO.parse(filename, "fasta")
    for rec in records:
        if len(rec.seq) > max_len:
            max_len = len(rec.seq)
            out = rec

    return out


def get_seq_lens(filename, seq_type):
    lens = []
    records = SeqIO.parse(filename, "fasta")
    for rec in records:
        alphabet = IUPAC.protein.letters if seq_type == 'prot' else IUPAC.unambiguous_dna.letters
        rec.seq = trim_sequence(rec.seq, alphabet)
        lens.append(len(rec.seq))

    return lens


def trim_sequence(seq, alphabet):
    end_pos = len(seq)
    for i in range(len(seq)):
        if seq[i] not in alphabet:
            end_pos = i
            break

    return seq[0:end_pos]


def save_to_fasta_file(records, output_file, mode="w+"):
    with open(output_file, mode) as output_file:
        SeqIO.write(records, output_file, "fasta")


def prepare_flank_sequences(seq_records, flank, ids=None):

    seq_recs = []
    seq_ranges = []
    seq_original_ranges = []
    for i, rec in enumerate(seq_records):
        seq_len = len(rec)
        seq_rec = rec.get_sequence(flank=flank)
        seq_original_range = rec.get_flank_range(flank=flank)
        flank_lens = rec.get_flank_lengths(flank)
        if ids:
            seq_rec.id = seq_rec.id + "_" + ids[i]

        seq_range = (flank_lens[0] + 1, flank_lens[0] + seq_len)
        seq_recs.append(seq_rec)
        seq_ranges.append(seq_range)
        seq_original_ranges.append(seq_original_range)

    return seq_recs, seq_ranges, seq_original_ranges
