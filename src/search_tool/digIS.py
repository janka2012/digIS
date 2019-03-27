from copy import copy
from copy import deepcopy

from ..blast.Blast import Blast
from ..blast.BlastX import BlastX
from ..blast.BlastN import BlastN
from ..common.csv_utils import write_csv
from ..common.sequence import *
from ..common.genbank import read_gb
from .RecordDigIS import RecordDigIS
from .digISClassifier import digISClassifier
from ..genbank.RecordGenbank import RecordGenbank
from ..hmmer.Hmmer import Hmmer


class digIS:

    def __init__(self, sequence_file, output_dir, config, genbank_file=None):
        self.sequence = sequence_file
        self.genbank_file = genbank_file
        self.sequence_len = get_seqlen(self.sequence)
        self.config = config
        self.genome = os.path.splitext(os.path.basename(sequence_file))[0]
        self.database = os.path.splitext(self.sequence)[0] + ".pep"
        translate_dna_seq_biopython(sequence=self.sequence, outseq=self.database)
        self.hmm = Hmmer(config.models, self.database)
        self.recs = []
        self.merged_records_idx = []
        self.output_dir = output_dir
        self.hmmer_output = os.path.join(output_dir, self.genome + '.hmmer3')
        self.csv_output = os.path.join(output_dir, self.genome + '.csv')
        self.filter_log = []
        self.genbank_overlap = []
        self.matched_recs = []
        self.fp = 0
        self.fn = 0

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

    def merge(self):
        """ Merging hits in particular distance """
        records_indexes = set(range(len(self.recs)))

        print("Number of records before merging: {}.".format(len(self.recs)))
        merged_records_indexes = self.merge_records(records_indexes)
        merged_records = [self.recs[rec_index] for rec_index in merged_records_indexes]
        self.recs = merged_records
        print("Number of records after merging: {}.".format(len(self.recs)))

    def merge_records(self, records_indexes):

        merged_records_indexes = self.__merge_records_in_distance(records_indexes)
        merged_records_indexes = self.__merge_overlapping_records(merged_records_indexes)

        return merged_records_indexes

    def __merge_records_in_distance(self, records_indexes):

        records_indexes_dc = deepcopy(records_indexes)
        for current_record_idx in copy(records_indexes_dc):
            if current_record_idx in records_indexes_dc:
                records_indexes_without_current_record = copy(records_indexes_dc)
                records_indexes_without_current_record.discard(current_record_idx)

                for other_record_idx in records_indexes_without_current_record:
                    current_record = self.recs[current_record_idx]
                    other_record = self.recs[other_record_idx]

                    if current_record.should_be_merged(other_record, self.config.max_merge_distance):
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

    def classification(self):

        # Extend hits with flank regions and get best Blast hits in original range
        classifier = digISClassifier()
        if self.recs:
            orf_blast_hits = Blast.get_best_blast_hits_in_range(self.recs, BlastX,
                                                                self.config.context_size_orf,
                                                                self.config.isfinder_orf_db)
            is_dna_blast_hits = Blast.get_best_blast_hits_in_range(self.recs, BlastN,
                                                                   self.config.context_size_is,
                                                                   self.config.isfinder_is_db)

            if self.genbank_file:
                genbank_recs = list(RecordGenbank(i, self.genome, "chr", self.sequence, self.sequence_len)
                                    for i in read_gb(self.genbank_file))
                genbank_recs = list(rec for rec in genbank_recs if rec.type not in ['source'])
                self.__assigns_genbank_annotation_by_overlap(genbank_recs, ignore_strand=True)
                ds_genbank_recs = self.__get_digis_genbank_map(genbank_recs)

                for digis_hit, gb_rec, orf_blast_hit, is_blast_hit in zip(self.recs, ds_genbank_recs,
                                                                          orf_blast_hits, is_dna_blast_hits):

                    classifier.classify(digis_hit, gb_rec, orf_blast_hit, is_blast_hit)
            else:
                for digis_hit, orf_blast_hit, is_blast_hit in zip(self.recs, orf_blast_hits, is_dna_blast_hits):
                    classifier.classify(digis_hit, None, orf_blast_hit, is_blast_hit)

    def __assigns_genbank_annotation_by_overlap(self, genbank_records=None, ignore_strand=False):

        unbinded_genbank_idx = set(list(range(len(genbank_records))))
        unbinded_digis_idx = set(list(range(len(self.recs))))
        match_idx = []

        # For each genbank record
        for gb_rec_idx, gb_rec in enumerate(genbank_records):
            # For each digis hit record
            max_overlap = 0
            for digis_rec_idx, digis_rec in enumerate(self.recs):
                # Test for overlap

                overlap = gb_rec.get_overlap_length(digis_rec, ignore_strand)
                if overlap > max_overlap:
                    max_overlap = overlap
                if overlap >= self.config.min_gb_overlap:
                    match_idx.append((gb_rec_idx, digis_rec_idx))
                    unbinded_genbank_idx.discard(gb_rec_idx)
                    unbinded_digis_idx.discard(digis_rec_idx)
            self.genbank_overlap.append(max_overlap)

        self.matched_recs = match_idx
        self.fn = len(unbinded_genbank_idx)
        self.fp = len(unbinded_digis_idx)

    def __get_digis_genbank_map(self, gb_recs=None):
        tool_recs = [[]] * len(self.recs)
        digis_match_idx = set(digis_idx for genbank_idx, digis_idx in self.matched_recs)
        for idx in digis_match_idx:
            for genbank_idx, digis_idx in self.matched_recs:
                if digis_idx == idx:
                    tool_recs[digis_idx] = tool_recs[digis_idx] + [gb_recs[genbank_idx]]
        return tool_recs

    def export(self, filename=None):
        csv_row = []
        csv_output = filename if filename else self.csv_output
        for rec in self.recs:
            csv_row.append([rec.qid, rec.sid, rec.qstart, rec.qend, rec.start, rec.end, rec.strand, rec.acc])
        write_csv(csv_row, csv_output, ["qid", "sid", "qstart", "qend", "sstart", "send", "strand", "acc"])

    def run(self, search=True, classification=True, export=True, debug=False):
        if search:
            self.search()
        self.parse()
        if debug:
            self.export(os.path.join(self.output_dir, self.genome + '_nonfilter.csv'))
        self.merge()
        if classification:
            self.classification()
        if export:
            self.export()

    def __str__(self):
        return '\n'.join(list(str(rec) for rec in self.recs))
