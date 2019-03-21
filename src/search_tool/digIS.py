import tempfile

from copy import copy
from copy import deepcopy


from src.blast.BlastX import BlastX
from src.common.csv_utils import write_csv
from src.common.misc import get_best_blast_hits_in_range, delete_file
from src.common.sequence import *
from src.search_tool.RecordDigIS import RecordDigIS
from src.hmmer.Hmmer import Hmmer

from definitions import ROOT_DIR


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

    def filter(self):
        """
        Filter merged records based on similarity with
        - positive sequence database (ISFinder ORF)
        - negative sequence database ()
        - negative profile database ()
        """

        self.__filter_by_sequence_db()
        # self.__filter_by_negative_profile_db()

        print(self.genome, 'Filtered', len(self.hmm.hsps) - len(self.recs))

        # Write filter log file
        with open(os.path.join(self.output_dir, self.genome + "_filter.log"), "w") as f:
            f.write("\n".join(self.filter_log))

    def __filter_by_sequence_db(self):

        if self.recs:

            positive_db_blast_hits = get_best_blast_hits_in_range(self.recs, BlastX,
                                                                  self.config.context_size_orf,
                                                                  self.config.isfinder_orf)
            negative_db_blast_hits = get_best_blast_hits_in_range(self.recs, BlastX,
                                                                  self.config.context_size_orf,
                                                                  self.config.negative_orf)

            records_indexes = set(range(len(self.recs)))
            print("Number of merged records before sequence db filtering: {}.".format(len(self.recs)))
            for i, record_idx in enumerate(copy(records_indexes)):
                if positive_db_blast_hits[i].score < negative_db_blast_hits[i].score:
                    print("Removing idx: {}.".format(record_idx))
                    records_indexes.discard(record_idx)
                    self.filter_log.append("Similar with negative database: " + str(negative_db_blast_hits[i]))

            # Update records after filtration
            save_recs = [self.recs[rec_index] for rec_index in records_indexes]
            self.recs = save_recs

        print("Number of merged records after sequence db filtering: {}.".format(len(self.recs)))

    def __filter_by_negative_profile_db(self):

        if self.recs:

            negative_profiles_db = os.path.join(ROOT_DIR, "data", "Hmmer_db", "negative_profiles_db.hmm")
            genome_id = os.path.splitext(os.path.basename(self.sequence))[0]
            hmmer_output = os.path.join(ROOT_DIR, "data", "Hmmer_db", genome_id + "_negative_profiles_db.hmmer")

            _, flanked_seqs_tmp_file = tempfile.mkstemp(prefix="digIS_", suffix=".fasta")

            ids = ["hit" + str(i) for i in range(len(self.recs))]

            flanked_seqs, flanked_seqs_ranges = prepare_flank_sequences(self.recs, self.config.context_size_orf, ids=ids)
            save_to_fasta_file(flanked_seqs, flanked_seqs_tmp_file)

            _, flanked_translated_seq_tmp_file = tempfile.mkstemp(prefix="digIS_", suffix=".pep")
            translate_dna_seq_biopython(sequence=flanked_seqs_tmp_file, outseq=flanked_translated_seq_tmp_file)
            # translate_dna_seq_biopython(sequence=flanked_seqs_tmp_file, outseq=flanked_translated_seq_tmp_file)

            evalue_threshold = 0.01
            hmmer = Hmmer(hmm=negative_profiles_db, seqfile=flanked_translated_seq_tmp_file)
            hmmer.run(tool="hmmscan", outfile=hmmer_output, evalue=evalue_threshold)
            hmmer.parse()

            print("Number of hmmscan hits: {}".format(len(hmmer.hits)))
            print("Number of false positive records before negative profile db filtering: {}".format(len(self.recs)))

            records_indexes = set(range(len(self.recs)))
            for hit in hmmer.hits:
                hit_id = int(hit.query_id.strip().split("_")[-2].replace("hit", ""))
                hit_dna_range = flanked_seqs_ranges[hit_id]
                best_hit = hit.get_best_hsp_in_range(self.sequence_len, hit_dna_range)

                if best_hit:
                    records_indexes.discard(hit_id)

            # Update records after filtration
            save_recs = [self.recs[rec_index] for rec_index in records_indexes]
            self.recs = save_recs

            print("Number of false positive records after negative profile db filtering: {}".format(len(self.recs)))

            # tempfiles clean up
            delete_file(flanked_seqs_tmp_file)
            delete_file(flanked_translated_seq_tmp_file)

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
        self.filter()
        if export:
            self.export()

    def __str__(self):
        return '\n'.join(list(str(rec) for rec in self.recs))
