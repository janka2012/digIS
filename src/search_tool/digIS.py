import shutil
from copy import copy
from copy import deepcopy
from math import ceil, log10

from ..common.gff_utils import write_gff
from ..common.classification import classification
from ..common.csv_utils import write_csv
from ..common.sequence import *
from ..hmmer.Hmmer import Hmmer
from ..genbank.RecordGenbank import RecordGenbank
from .RecordDigIS import RecordDigIS


class digIS:

    def __init__(self, config, genome, genbank_features):
        self.config = config
        self.genome = genome
        self.genbank_features = genbank_features
        self.hmmer = Hmmer()
        self.hmmsearch_output = os.path.join(self.config.output_dir, "hmmer", str(self.genome.name) + '_hmmsearch.hmmer3')
        self.phmmer_output = os.path.join(self.config.output_dir, "hmmer", str(self.genome.name) + '_phmmer.hmmer3')
        self.log_output = os.path.join(self.config.output_dir, "logs", str(self.genome.name) + '_filter.log')
        self.output = os.path.join(self.config.output_dir, "results", str(self.genome.name) + "." + str(self.config.out_format))
        self.recs = []
        self.genbank_overlap = []
        self.matched_recs = []
        self.filter_log = []
        self.classifier_recs = []

    def search(self):
        print("Hmmer search...")
        self.search_models()
        self.search_outliers()

    def search_models(self):
        self.hmmer.run(tool="hmmsearch", hmmfile=self.config.models, seqdb=self.genome.orf_db,
                       outfile=self.hmmsearch_output, curated_models=True)

    def search_outliers(self):
        print("Searching using Hmmer...")
        self.hmmer.run(tool="phmmer", hmmfile=self.config.outliers_fasta, seqdb=self.genome.orf_db,
                       outfile=self.phmmer_output, evalue="0.001", cevalue="0.001")

    def parse(self):
        print("Parsing Hmmer outputs...")
        self.parse_hmmer_output(self.hmmsearch_output)
        self.parse_hmmer_output(self.phmmer_output)

    def parse_hmmer_output(self, hmmer_output):
        new_recs_added = self.hmmer.parse(hmmer_output)

        if new_recs_added:
            for hsp in self.hmmer.hsps:
                hit_range = (hsp.sstart, hsp.send)
                frame = int(hsp.sid.strip()[-1])
                sid = hsp.sid[:-2]
                strand = "+" if frame <= 3 else "-"
                dna_range = transform_range(hit_range[0], hit_range[1], frame, self.genome.length)
                self.recs.append(RecordDigIS.from_hmmer(hsp, sid, dna_range[0], dna_range[1], strand, self.genome.name,
                                                        "chr", self.genome.seq, self.genome.length))

    def merge(self):
        """ Merging hits in particular distance """
        records_indexes = set(range(len(self.recs)))

        print("Merging hits...")
        print("Number of records before merging: {}.".format(len(self.recs)))
        merged_records_indexes = self.merge_records(records_indexes)
        merged_records = [self.recs[rec_index] for rec_index in merged_records_indexes]
        self.recs = merged_records
        print("Number of records after merging: {}.".format(len(self.recs)))

        # Write filter log file
        with open(self.log_output, "w+") as f:
            f.write("\n".join(self.filter_log))

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

    def filter(self):
        print("Filtering hits shorter or equal than {} bp".format(self.config.min_hit_length))
        recs_filt = [rec for rec in self.recs if abs(rec.start - rec.end + 1) >= self.config.min_hit_length]
        self.recs = recs_filt

    def classification(self):
        if self.config.genbank_file:
            print("Classification with GenBank annotation...")
            genbank_recs = list(RecordGenbank(i, self.genome.name, "chr", self.genome.seq, self.genome.length)
                                for i in self.genbank_features)
            genbank_recs = list(rec for rec in genbank_recs if rec.type not in ['source'])
        else:
            print("Classification without GenBank annotation...")
            genbank_recs = None

        self.classifier_recs = classification(self.recs, genbank_recs, self.config.context_size_orf,
                                              self.config.context_size_is, self.config.min_gb_overlap,
                                              self.config.isfinder_orf_db, self.config.isfinder_is_db)

    def export_records(self):
        csv_row = []
        csv_header = ["id", "level", "qid", "sid", "qstart", "qend", "sstart", "send", "strand", "acc"]
        id_width = ceil(log10(len(self.recs))) if len(self.recs) > 0 else 1

        for i, rec in enumerate(self.recs):
            if self.classifier_recs[i].similarity_is in ['medium', 'strong']:
                rec.start = self.classifier_recs[i].blast_is_dna.query_start
                rec.end = self.classifier_recs[i].blast_is_dna.query_end
                level = 'is'
            elif self.classifier_recs[i].similarity_orf in ['medium', 'strong']:
                rec.start = self.classifier_recs[i].blast_orf.query_start
                rec.end = self.classifier_recs[i].blast_orf.query_end
                level = 'orf'
            else:
                level = 'domain'

            id = '_'.join([rec.sid, str(i).zfill(id_width), level])
            search_header, search_row = rec.to_csv()
            class_header, class_row = self.classifier_recs[i].to_csv()

            header = ['id', 'level'] + search_header + class_header
            row = [id, level] + search_row + class_row

            csv_row.append(row)
            csv_header = header
        return csv_header, csv_row

    def export_summary_stats(self):
        # Create Family dictionary
        family_dict = {}
        for rec in self.recs:
            id = 'others' if '_' in rec.qid else rec.qid
            if id in family_dict.keys():
                num, bps = family_dict[id]
                family_dict[id] = (num+1, bps+len(rec))
            else:
                family_dict[id] = (1, len(rec))

        # Create list of summary records
        sum_recs = []
        genome_len = self.genome.length
        total_num = 0
        total_bps = 0
        total_pdna = 0.0
        fam_id_sorted = list(family_dict.keys())
        fam_id_sorted.sort()
        for fam_id in fam_id_sorted:
            fam_num, fam_bps = family_dict[fam_id]
            fam_pdna = fam_bps / genome_len * 100
            sum_recs.append([self.genome.name, fam_id, fam_num, fam_bps, genome_len, fam_pdna])
            total_num += fam_num
            total_bps += fam_bps
            total_pdna += fam_pdna
        sum_recs.append([self.genome.name, 'total', total_num, total_bps, genome_len, total_pdna])

        return sum_recs

    def export(self, filename=None):
        print("Exporting output...")
        output = filename if filename else self.output
        csv_header, csv_row = self.export_records()
        if self.config.out_format == "csv":
            write_csv(csv_row, output, csv_header)
        elif self.config.out_format == "gff":
            write_gff(csv_row, output, csv_header)

    def run(self, search=True, export=False):
        print("===== Processing of", self.genome.name, "sequence =====")
        if search:
            self.search()
        self.parse()
        self.merge()
        self.filter()
        self.classification()
        if export:
            self.export()

    def __str__(self):
        return '\n'.join(list(str(rec) for rec in self.recs))
