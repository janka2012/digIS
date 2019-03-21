import subprocess
import sys

from Bio import SearchIO, SeqIO

from src.common.csv_utils import write_csv
from src.common.misc import check_if_file_exists, change_path_to_linux
from src.common.misc import check_evalue
from src.common.sequence import get_sequence_record_ids
from src.hmmer.HmmerHit import HmmerHit


class Hmmer:
    def __init__(self, hmm, seqfile):
        """ Create a new Hmmer instance.

        model       HMM model used in hmmsearch (e.g. "IS1_cut.hmm")
        database    fasta file, used as database in hmmsearch (e.g. "test.fasta")
        hits        list of hmmer hits
        hsps        list of hmmer hsps
        outfile     file, where output from hmmsearch is stored (e.g. "out.hmmer")
        """
        self.hmm = check_if_file_exists(hmm)
        self.seqfile = check_if_file_exists(seqfile)
        self.hits = []
        self.hsps = []
        self.outfile = ""

    def run(self, tool, outfile, evalue="0.001", curated_models=False):
        check_evalue(evalue)

        if outfile:
            self.outfile = outfile
            cmd = self.__build_command(tool=tool, evalue=str(evalue), outfile=outfile, curated_models=curated_models)
            if sys.platform == 'win32':
                cmd = ['bash.exe', '-c', ' '.join(cmd)]
            self.__run_tool(cmd)
        else:
            raise AttributeError("Output file argument is required.")

    def parse(self, outfile=None):

        hmmer_file = None
        try:
            hmmer_file = outfile if outfile else self.outfile
            check_if_file_exists(hmmer_file)
        except FileNotFoundError:
            print("No hmmer output file set.")

        self.hits = []
        self.hsps = []
        for res in SearchIO.parse(hmmer_file, 'hmmsearch3-domtab'):
            for hit in res.hits:
                print(hit)
                self.hits.append(HmmerHit(hit, res.seq_len))
        for hit in self.hits:
            self.hsps += hit.hsps

    def save_hmmer_output_to_csv(self, output_csv, hmmer_outfile=None):
        hmmer_file = None
        try:
            hmmer_file = hmmer_outfile if hmmer_outfile else self.outfile
            check_if_file_exists(hmmer_file)
        except FileNotFoundError:
            print("No hmmer output file set.")

        header = ["subject_id", "subject_len", "query_id", "query_len", "seq_evalue", "seq_bitscore", "seq_bias",
                  "dom_idx", "dom_num", "dom_cevalue", "dom_evalue", "dom_bitscore", "dom_bias", "query_start",
                  "query_end", "subject_start", "subject_end", "subject_env_start", "subject_env_end", "acc_avg",
                  "subject_description"]
        rows = []
        for res in SearchIO.parse(hmmer_file, 'hmmsearch3-domtab'):
            for hit in res.hits:
                for hsp in hit.hsps:

                    row = [hit.id, hit.seq_len, res.id, res.seq_len, hit.evalue, hit.bitscore, hit.bias,
                           hsp.domain_index, len(hit.hsps), hsp.evalue_cond, hsp.evalue, hsp.bitscore, hsp.bias,
                           hsp.query_start+1, hsp.query_end, hsp.hit_start+1, hsp.hit_end, hsp.env_start,
                           hsp.env_end, hsp.acc_avg, hit.description]

                    rows.append(row)
        write_csv(rows, output_csv, header)

    def save_hmmer_hits_to_fasta(self, output):
        import os
        ids = [hit.subject_id for hit in self.hits]
        print(os.path.exists(self.seqfile))
        recs = get_sequence_record_ids(self.seqfile, ids)
        print("Writing")
        SeqIO.write(recs, output, "fasta")

    def str_hsps(self):
        return '\n'.join(list(str(i) for i in self.hsps))

    def __build_command(self, tool, evalue, outfile, curated_models):
        """
        tool [options] <hmmdb> <seqfile>
        """

        cmd = [tool, "--noali"]

        if curated_models:
            cmd.extend(["--cut_nc"])
        else:
            cmd.extend(["-E", evalue])

        if sys.platform == 'win32':
            outfile = change_path_to_linux(outfile)
            hmm = change_path_to_linux(self.hmm)
            seqfile = change_path_to_linux(self.seqfile)
            cmd.extend(["--domtblout", outfile, hmm, seqfile])
        else:
            cmd.extend(["--domtblout", outfile, self.hmm, self.seqfile])
        return cmd

    @staticmethod
    def __run_tool(cmd):

        try:
            subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)
        except subprocess.CalledProcessError as e:
            print("An error occurred when calling {}.".format(cmd))
            print(e)

    def __len__(self):
        return len(self.hits)

    def __str__(self):
        return '\n'.join(list(str(i) for i in self.hits))
