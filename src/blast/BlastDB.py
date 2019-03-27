import os

from subprocess import call
from Bio import SeqIO

from definitions import ROOT_DIR


class BlastDB:
    def __init__(self, fasta_file, db_name, dir_path="", db_type="prot"):
        self.fasta_file = fasta_file
        self.dir_path = dir_path
        self.db_name = db_name.replace(" ", "_")
        self.db_type = db_type
        self.db_output_dir = os.path.join(ROOT_DIR, "data", "Blast_db")
        self.db_output_file = os.path.join(self.db_output_dir, self.db_name)

    def make_db(self):
        if os.path.exists(self.fasta_file):
            call(["makeblastdb", "-in", self.fasta_file, "-parse_seqids",
                  "-dbtype", self.db_type, "-out", self.db_output_file])
        else:
            raise FileExistsError("File {} does not exists.".format(self.fasta_file))

    def prepare_fasta_file(self, out_fasta):
        all_records = []

        for fasta_file in self.__get_fasta_files():
            collected_records = self.__collect_fasta_records_from_file(fasta_file)
            all_records.extend(collected_records)
            SeqIO.write(all_records, out_fasta, "fasta")

    def __get_fasta_files(self):
        if os.path.exists(self.dir_path) and os.path.isdir(self.dir_path):
            for (dirpath, dirnames, filenames) in os.walk(self.dir_path):
                for filename in filenames:
                    file = os.path.join(dirpath, filename)
                    if os.path.exists(file) and file.endswith(".fasta"):
                        yield file

    def __collect_fasta_records_from_file(self, file):
        formatted_records = []
        for rec in SeqIO.parse(file, "fasta"):
            self.__format_record_header(rec)
            formatted_records.append(rec)
        return formatted_records

    def __format_record_header(self, record):
        record.id = "gnl|" + self.db_name + "|" + record.description.replace("|", ";")
        record.description = ""
