import shutil
from copy import copy
from copy import deepcopy
from math import ceil, log10

from ..common.gff_utils import write_gff
from ..common.classification import classification
from ..common.Classifier import Classifier
from ..common.csv_utils import write_csv
from ..common.sequence import *
from ..hmmer.Hmmer import Hmmer
from ..genbank.RecordGenbank import RecordGenbank
from .RecordDigIS import RecordDigIS
from .RecordDigISAttrib import RecordDigISAttrib
from ..blast.Blast import Blast
from ..blast.BlastX import BlastX
from ..blast.BlastN import BlastN


class digIS:

    def __init__(self, config, genome, genbank_features):
        self.config = config
        self.genome = genome
        self.genbank_features = genbank_features
        self.hmmer = Hmmer()
        self.hmmsearch_output = os.path.join(self.config.output_dir, "hmmer", str(self.genome.name) + '_hmmsearch.hmmer3')
        self.annotated_output = os.path.join(self.config.output_dir, "hmmer", str(self.genome.name) + '_annotated.csv')
        self.phmmer_output = os.path.join(self.config.output_dir, "hmmer", str(self.genome.name) + '_phmmer.hmmer3')
        self.log_output = os.path.join(self.config.output_dir, "logs", str(self.genome.name) + '_filter.log')
        self.output_csv = os.path.join(self.config.output_dir, "results", str(self.genome.name) + ".csv")
        self.output_gff = os.path.join(self.config.output_dir, "results", str(self.genome.name) + ".gff")
        self.recs = []
        self.recs_attrib = []
        # self.extension_level = []
        # self.genbank_overlap = []
        # self.matched_recs = []
        self.filter_log = []
        # self.classifier_recs = []

    def __get_valid_indexes_and_records(self):
        indexes = []
        records = []
        for i, rec_attrib in enumerate(self.recs_attrib):
            if rec_attrib.status == 'valid':
                indexes.append(i)
                records.append(self.recs[i])
        return indexes, records

    def search(self):
        print("Seed search...")
        self.search_models()
        self.search_outliers()

    def search_models(self):
        self.hmmer.run(tool="hmmsearch", hmmfile=self.config.models, seqdb=self.genome.orf_db,
                       outfile=self.hmmsearch_output)

    def search_outliers(self):
        self.hmmer.run(tool="phmmer", hmmfile=self.config.outliers_fasta, seqdb=self.genome.orf_db,
                       outfile=self.phmmer_output)

    def parse(self):
        print("Parsing Hmmer outputs...")
        self.parse_hmmer_output('model', self.hmmsearch_output)
        self.parse_hmmer_output('outlier', self.phmmer_output)

    def parse_hmmer_output(self, source_type, hmmer_output):
        new_recs = self.hmmer.parse(hmmer_output)

        if new_recs:
            for hsp in self.hmmer.flat_hsps:
                hit_range = (hsp.sstart, hsp.send)
                frame = int(hsp.sid.strip()[-1])
                sid = hsp.sid[:-2]
                strand = "+" if frame <= 3 else "-"
                dna_range = transform_range(hit_range[0], hit_range[1], frame, self.genome.length)
                self.recs.append(RecordDigIS.from_hmmer(hsp, sid, dna_range[0], dna_range[1], strand, self.genome.name,
                                                        "chr", self.genome.seq, self.genome.length))
                self.recs_attrib.append(RecordDigISAttrib(source_type, hsp))

    def merge(self):
        """ Merging hits in particular distance """
        records_indexes, _records = self.__get_valid_indexes_and_records()

        print("Seed Merging...")
        print("Number of records before merging: {}.".format(len(records_indexes)))
        merged_records_indexes = self.merge_records(set(records_indexes))
        print("Number of records after merging: {}.".format(len(merged_records_indexes)))

    def merge_records(self, records_indexes):

        merged_records_indexes = self.__merge_records_in_distance(records_indexes)
        merged_records_indexes = self.__merge_overlapping_records(merged_records_indexes)

        return merged_records_indexes

    def __merge_records_in_distance(self, records_indexes):

        records_indexes_copy = records_indexes.copy()
        for current_record_idx in records_indexes:
            if current_record_idx in records_indexes_copy:
                records_indexes_without_current_record = records_indexes_copy.copy()
                records_indexes_without_current_record.discard(current_record_idx)

                for other_record_idx in records_indexes_without_current_record:
                    current_record = self.recs[current_record_idx]
                    other_record = self.recs[other_record_idx]
                    other_record_attrib = self.recs_attrib[other_record_idx]

                    if current_record.should_be_merged(other_record, self.config.max_merge_distance):
                        spol = current_record.start < other_record.start
                        qpol = current_record.qstart < other_record.qstart

                        if (current_record.strand == '+' and spol == qpol) or (current_record.strand == '-' and spol != qpol):
                            self.filter_log.append("{}:{}  {} merged with neighbouring element: {}".format(current_record_idx, other_record_idx, current_record, other_record))
                            current_record.merge(other_record, merge_type="distance")
                            records_indexes_copy.discard(other_record_idx)
                            other_record_attrib.status = 'merged_distance'

        return records_indexes_copy

    def __merge_overlapping_records(self, records_indexes):

        records_indexes_copy = records_indexes.copy()

        for current_record_idx in records_indexes:
            if current_record_idx in records_indexes_copy:
                records_indexes_without_current_record = records_indexes_copy.copy()
                records_indexes_without_current_record.discard(current_record_idx)

                for other_record_idx in records_indexes_without_current_record:
                    current_record = self.recs[current_record_idx]
                    other_record = self.recs[other_record_idx]
                    other_record_attrib = self.recs_attrib[other_record_idx]

                    if current_record.get_overlap_length(other_record) > 0:
                        self.filter_log.append("{}:{}  {} merged with overlapping element: {}".format(current_record_idx, other_record_idx, current_record, other_record))
                        current_record.merge(other_record, merge_type="overlap")
                        records_indexes_copy.discard(other_record_idx)
                        other_record_attrib.status = 'merged_overlap'

        return records_indexes_copy

    def seed_extension(self):
        print("Seed Extension...")
        indexes, records = self.__get_valid_indexes_and_records()
        if len(indexes) > 0:
            orf_blast_pos = Blast.get_max_blast_hits_in_range(records, BlastX, self.config.context_size_orf, self.config.isfinder_orf_db,
                                                              min_overlap=1, positive_subject_strand_only=True)
            is_blast_pos = Blast.get_max_blast_hits_in_range(records, BlastN, self.config.context_size_is, self.config.isfinder_is_db,
                                                             min_overlap=1, positive_subject_strand_only=True)

            for idx, rec, orf_pos, is_pos in zip(indexes, records, orf_blast_pos, is_blast_pos):
                # Extension at the level of ORF
                orf_start, orf_end = orf_pos
                if not (orf_start == 0 and orf_end == 0):
                    if orf_start < rec.start or orf_end > rec.end:
                        rec.start = min(rec.start, orf_start)
                        rec.end = max(rec.end, orf_end)
                        self.recs_attrib[idx].extension_level = 'orf'
                # Extension at the level of IS
                is_start, is_end = is_pos
                if not (is_start == 0 and is_end == 0):
                    if is_start < rec.start or is_end > rec.end:
                        rec.start = min(rec.start, is_start)
                        rec.end = max(rec.end, is_end)
                        self.recs_attrib[idx].extension_level = 'is'

    def filter_by_length(self):
        print("Filtering hits shorter or equal than {} bp".format(self.config.min_hit_length))
        for idx, (rec_attrib, rec) in enumerate(zip(self.recs_attrib, self.recs)):
            if rec_attrib.status == 'valid' and len(rec) < self.config.min_hit_length:
                rec_attrib.status = 'filtered_length'
                self.filter_log.append("{}: {} filtered because of the length".format(idx, rec))

    def filter_by_evalue(self, evalue):
        print("Filtering hits with e-value higher than {}.".format(evalue))
        for idx, (rec_attrib, rec) in enumerate(zip(self.recs_attrib, self.recs)):
            if rec_attrib.status == 'valid' and rec.evalue > evalue:
                rec_attrib.status = 'filtered_evalue'
                self.filter_log.append("{}: {} filtered because of the low e-value".format(idx, rec))


    def filter_by_score(self):
        print("Filtering hits by cut off thresholds.")
        cutoffs_dict = self.__load_threshold()

        for idx, (rec_attrib, rec) in enumerate(zip(self.recs_attrib, self.recs)):
            model = rec.qid
            if model in cutoffs_dict.keys():
                score = cutoffs_dict[model]
                if rec.score < score:
                    rec_attrib.status = 'filtered_score'
                    self.filter_log.append("{}: {} filtered because of the low score".format(idx, rec))

    def __load_threshold(self):
        import os
        import csv
        from definitions import ROOT_DIR
        nc_thresholds = os.path.join(ROOT_DIR, "digIS", "data", "models", "hmm", "thresholds.txt")

        with open(nc_thresholds) as tsvfile:
            cutoffs_dict = {}
            reader = csv.reader(tsvfile, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
            next(reader)
            for row in reader:
                cutoffs_dict[str(row[0])] = float(row[2])

        return cutoffs_dict

    def classification(self):
        if self.config.genbank_file:
            print("Classification with GenBank annotation...")
            genbank_recs = list(RecordGenbank(i, self.genome.name, "chr", self.genome.seq, self.genome.length)
                                for i in self.genbank_features)
            genbank_recs = list(rec for rec in genbank_recs if rec.type not in ['source'])
        else:
            print("Classification without GenBank annotation...")
            genbank_recs = None

        indexes, records = self.__get_valid_indexes_and_records()
        classifier_recs = classification(records, genbank_recs, self.config.context_size_orf,
                                              self.config.context_size_is, self.config.min_gb_overlap,
                                              self.config.isfinder_orf_db, self.config.isfinder_is_db)

        for idx, class_rec in zip(indexes, classifier_recs):
            self.recs_attrib[idx].classification = class_rec

    def export_records(self):
        csv_row = []
        csv_header = ["id", "level", "qid", "sid", "qstart", "qend", "sstart", "send", "strand", "acc"]
        id_width = ceil(log10(len(self.recs))) if len(self.recs) > 0 else 1
        indexes, records = self.__get_valid_indexes_and_records()

        search_header = RecordDigIS.get_csv_header()
        class_header = Classifier.get_csv_header()
        csv_header = ['id', 'level'] + search_header + class_header

        for i, (idx, rec) in enumerate(zip(indexes, records)):
            extension_level = self.recs_attrib[idx].extension_level
            class_rec = self.recs_attrib[idx].classification

            id = '_'.join([rec.sid, str(i).zfill(id_width), extension_level])
            row = [id, extension_level] + rec.to_csv() + class_rec.to_csv()
            csv_row.append(row)

        return csv_header, csv_row

    # For debug purposess only
    def export_annotated_records(self):
        csv_header = RecordDigISAttrib.get_csv_header() + RecordDigIS.get_csv_header()
        csv_rows = []
        for rec_attrib, rec in zip(self.recs_attrib, self.recs):
            csv_rows.append(rec_attrib.to_csv() + rec.to_csv())

        write_csv(csv_rows, self.annotated_output, csv_header)


    def export_summary_stats(self):
        _indexes, records = self.__get_valid_indexes_and_records()
        # Create Family dictionary
        family_dict = {}
        for rec in records:
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

    def export(self):
        print("Exporting output...")
        csv_header, csv_row = self.export_records()
        write_csv(csv_row, self.output_csv, csv_header)
        write_gff(csv_row, self.output_gff, csv_header)

    def run(self, search=True, export=True):
        print("===== Processing of", self.genome.name, "sequence =====")
        if search:
            self.search()
        self.parse()
        self.filter_by_evalue(self.config.outliers_evalue)
        self.merge()
        self.seed_extension()
        self.filter_by_score()
        self.filter_by_length()
        self.classification()
        if export:
            self.export()
            self.export_annotated_records()

        # Write output log file
        with open(self.log_output, "w+") as f:
            f.write("\n".join(self.filter_log))

    def __str__(self):
        return '\n'.join(list(str(rec) for rec in self.recs))
