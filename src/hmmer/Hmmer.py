import subprocess
import sys

from Bio import SearchIO, SeqIO

from ..common.csv_utils import write_csv
from ..common.misc import check_if_file_exists, change_path_to_linux
from ..common.sequence import get_sequence_record_ids
from ..hmmer.HmmerHit import HmmerHit


class Hmmer:
    def __init__(self):
        """ Create a new Hmmer instance.

        model       HMM model used in hmmsearch (e.g. "IS1_cut.hmm")
        database    fasta file, used as database in hmmsearch (e.g. "test.fasta")
        hits        list of hmmer hits
        hsps        list of hmmer hsps
        outfile     file, where output from hmmsearch is stored (e.g. "out.hmmer")
        """
        self.hmm = ""
        self.seqfile = ""
        self.hits = []
        self.hsps = []

    def run(self, tool, hmmfile, seqdb, outfile, evalue=None, cevalue=None, curated_models=False):

        if not check_if_file_exists(seqdb):
            print("File {} does not exist or is not a file.".format(seqdb))
            print("Try to run digIS with translate option turned on.")
            exit(1)

        if outfile:
            cmd = self.__build_command(tool=tool, hmmfile=hmmfile, seqdb=seqdb,
                                       outfile=outfile, curated_models=curated_models, evalue=evalue, cevalue=cevalue)
            if sys.platform == 'win32':
                cmd = ['bash.exe', '-c', ' '.join(cmd)]
            self.__run_tool(cmd)
        else:
            raise AttributeError("Output file argument is required.")

    def parse(self, outfile):

        new_recs_added = False

        try:
            check_if_file_exists(outfile)
        except FileNotFoundError:
            print("No hmmer output file set.")

        hmmer_res = list(SearchIO.parse(outfile, 'hmmsearch3-domtab'))
        if len(hmmer_res) > 0:
            new_recs_added = True
            for res in hmmer_res:
                for hit in res.hits:
                    self.hits.append(HmmerHit(hit, res.seq_len))

            for hit in self.hits:
                self.hsps += hit.hsps
        return new_recs_added

    def save_hmmer_output_to_csv(self, output_csv, hmmer_outfile):
        try:
            check_if_file_exists(hmmer_outfile)
        except FileNotFoundError:
            print("No hmmer output file set.")

        header = ["subject_id", "subject_len", "query_id", "query_len", "seq_evalue", "seq_bitscore", "seq_bias",
                  "dom_idx", "dom_num", "dom_cevalue", "dom_evalue", "dom_bitscore", "dom_bias", "query_start",
                  "query_end", "subject_start", "subject_end", "subject_env_start", "subject_env_end", "acc_avg",
                  "subject_description"]
        rows = []
        for res in SearchIO.parse(hmmer_outfile, 'hmmsearch3-domtab'):
            for hit in res.hits:
                for hsp in hit.hsps:

                    row = [hit.id, hit.seq_len, res.id, res.seq_len, hit.evalue, hit.bitscore, hit.bias,
                           hsp.domain_index, len(hit.hsps), hsp.evalue_cond, hsp.evalue, hsp.bitscore, hsp.bias,
                           hsp.query_start+1, hsp.query_end, hsp.hit_start+1, hsp.hit_end, hsp.env_start,
                           hsp.env_end, hsp.acc_avg, hit.description]

                    rows.append(row)
        write_csv(rows, output_csv, header)

    def save_hmmer_hits_to_fasta(self, output):
        ids = [hit.subject_id for hit in self.hits]
        recs = get_sequence_record_ids(self.seqfile, ids)
        SeqIO.write(recs, output, "fasta")

    def str_hsps(self):
        return '\n'.join(list(str(i) for i in self.hsps))

    def __build_command(self, tool, hmmfile, seqdb, outfile, curated_models, evalue=None, cevalue=None):
        """
        tool [options] <hmmdb> <seqfile>
        """

        cmd = [tool, "--noali"]

        if curated_models and (evalue or cevalue):
            raise ValueError("You can not set both - noise thresholds and evalue. Set either thresholds or evalue.")

        if evalue:
            cmd.extend(["-E", evalue])
            if cevalue:
                cmd.extend(["--domE", cevalue])

        elif curated_models:
            cmd.extend(["--cut_nc"])

        if sys.platform == 'win32':
            outfile = change_path_to_linux(outfile)
            hmmfile = change_path_to_linux(hmmfile)
            seqdb = change_path_to_linux(seqdb)

        cmd.extend(["--domtblout", outfile, hmmfile, seqdb])
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
