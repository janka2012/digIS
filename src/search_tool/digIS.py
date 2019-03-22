import tempfile

from copy import copy
from copy import deepcopy


from src.blast.BlastX import BlastX
from src.common.csv_utils import write_csv
from src.common.misc import get_best_blast_hits_in_range, delete_file
from src.common.sequence import *
from src.search_tool.RecordDigIS import RecordDigIS
from src.hmmer.Hmmer import Hmmer

from src.common.definitions import ROOT_DIR


class digIS:

    def __init__(self, model, sequence_file, output_dir, config):
        self.model = model
        self.sequence = sequence_file
        self.sequence_len = get_seqlen(self.sequence)
        self.config = config
        self.genome = os.path.splitext(os.path.basename(sequence_file))[0]
        self.database = os.path.splitext(self.sequence)[0] + ".pep"
        translate_dna_seq_biopython(sequence=self.sequence, outseq=self.database)
        self.hmm = Hmmer(self.model, self.database)
        self.recs = []
        self.merged_records_idx = []
        self.output_dir = output_dir
        self.hmmer_output = os.path.join(output_dir, self.genome + '.hmmer3')
        self.csv_output = os.path.join(output_dir, self.genome + '.csv')
        self.filter_log = []

    def search(self):
        self.hmm.run(tool="hmmsearch", outfile=self.hmmer_output, curated_models=True)

    def parse(self):
        self.hmm.parse(self.hmmer_output)
        self.recs = []
        for hsp in self.hmm.hsps:
            hit_range = (hsp.sstart, hsp.send)
            frame = int(hsp.sid.strip()[-1])
            sid = hsp.sid[:-2]
            strand = "+" if frame <= 3 else "-"
            dna_range = transform_range(hit_range[0], hit_range[1], frame, self.sequence_len)
            self.recs.append(RecordDigIS.from_hmmer(hsp, sid, dna_range[0], dna_range[1], strand,
                                                    self.genome, "chr", self.sequence, self.sequence_len))

    def merge(self, merge_distance=300):
        """ Merging hits in particular distance """
        records_indexes = set(range(len(self.recs)))

        print("Number of records before merging: {}.".format(len(self.recs)))
        merged_records_indexes = self.merge_records(records_indexes, merge_distance)
        merged_records = [self.recs[rec_index] for rec_index in merged_records_indexes]
        self.recs = merged_records
        print("Number of records after merging: {}.".format(len(self.recs)))
        # self.merged_records_idx = merged_records_indexes

    def merge_records(self, records_indexes, merge_distance=300):

        merged_records_indexes = self.__merge_records_in_distance(records_indexes, merge_distance)
        merged_records_indexes = self.__merge_overlapping_records(merged_records_indexes)

        return merged_records_indexes

    def __merge_records_in_distance(self, records_indexes, merge_distance):

        records_indexes_dc = deepcopy(records_indexes)
        for current_record_idx in copy(records_indexes_dc):
            if current_record_idx in records_indexes_dc:
                records_indexes_without_current_record = copy(records_indexes_dc)
                records_indexes_without_current_record.discard(current_record_idx)

                for other_record_idx in records_indexes_without_current_record:
                    current_record = self.recs[current_record_idx]
                    other_record = self.recs[other_record_idx]

                    if current_record.should_be_merged(other_record, merge_distance):
                        spol = current_record.start < other_record.start
                        qpol = current_record.qstart < other_record.qstart

                        if (current_record.strand == '+' and spol == qpol) or (current_record.strand == '-' and spol != qpol):
                            current_record.merge(other_record)
                            records_indexes_dc.discard(other_record_idx)
                            self.filter_log.append("{} merged with neighbouring element: {}".format(current_record, other_record))
        return records_indexes_dc

    def __merge_overlapping_records(self, records_indexes):

        records_indexes_dc = deepcopy(records_indexes)

        for current_record_idx in copy(records_indexes_dc):
            if current_record_idx in records_indexes_dc:
                records_indexes_without_current_record = copy(records_indexes_dc)
                records_indexes_without_current_record.discard(current_record_idx)
                for other_record_idx in records_indexes_without_current_record:
                    current_record = self.recs[current_record_idx]
                    other_record = self.recs[other_record_idx]
                    if current_record.get_overlap_length(other_record) > 0:
                        current_record.merge(other_record)
                        records_indexes_dc.discard(other_record_idx)
                        self.filter_log.append("{} merged with overlapping element: {}".format(current_record, other_record))
        return records_indexes_dc

    def export(self, filename=None):
        csv_row = []
        csv_output = filename if filename else self.csv_output
        for rec in self.recs:
            csv_row.append([rec.qid, rec.sid, rec.qstart, rec.qend, rec.start, rec.end, rec.strand, rec.acc])
        write_csv(csv_row, csv_output, ["qid", "sid", "qstart", "qend", "sstart", "send", "strand", "acc"])

    def run(self, search=True, export=True, debug=False):
        if search:
            self.search()
        self.parse()
        if debug:
            self.export(os.path.join(self.output_dir, self.genome + '_nonfilter.csv'))
        self.merge()
        if export:
            self.export()

    def __str__(self):
        return '\n'.join(list(str(rec) for rec in self.recs))
