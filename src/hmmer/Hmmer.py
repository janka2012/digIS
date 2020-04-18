import subprocess
import sys

from Bio import SearchIO, SeqIO

from definitions import NUM_THREADS
from ..common.csv_utils import write_csv
from ..common.misc import check_if_file_exists, change_path_to_linux
from ..common.sequence import get_sequence_record_ids
from ..hmmer.HmmerHspFlat import HmmerHspFlat


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
        self.flat_hsps = []

    def run(self, tool, hmmfile, seqdb, outfile, cevalue=None):

        if not check_if_file_exists(seqdb):
            print("File {} does not exist or is not a file.".format(seqdb))
            print("Try to run digIS with translate option turned on.")
            exit(1)

        if outfile:
            cmd = self.__build_command(tool=tool, hmmfile=hmmfile, seqdb=seqdb, outfile=outfile, cevalue=cevalue)
            if sys.platform == 'win32':
                cmd = ['bash.exe', '-c', ' '.join(cmd)]
            self.__run_tool(cmd)
        else:
            raise AttributeError("Output file argument is required.")

    def parse(self, outfile):
        self.flat_hsps = []
        new_recs = False

        try:
            check_if_file_exists(outfile)
        except FileNotFoundError:
            print("No hmmer output file set.")

        hmmer_res = list(SearchIO.parse(outfile, 'hmmsearch3-domtab'))
        if len(hmmer_res) > 0:
            new_recs = True
            for query in hmmer_res:
                for hit in query.hits:
                    for hsp in hit.hsps:
                        self.flat_hsps.append(HmmerHspFlat(query, hit, hsp))

        return new_recs

    def to_csv(self, output_csv=None):
        header = []
        rows = []
        for hsp in self.flat_hsps:
            header, row = hsp.to_csv()
            rows.append(row)
        
        if output_csv:
            write_csv(rows, output_csv, header)

        return header, rows

    def __build_command(self, tool, hmmfile, seqdb, outfile, cevalue=None):
        """
        tool [options] <hmmdb> <seqfile>
        """

        cmd = [tool, "--noali"]

        if NUM_THREADS != 0:
            cmd.extend(["--cpu", str(NUM_THREADS)])

        if cevalue:
            cmd.extend(["--domE", str(cevalue)])
        else:
            cmd.extend(["--domT", "0.0"])

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
