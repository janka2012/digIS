from Bio.Blast.Applications import NcbiblastnCommandline
from src.blast.Blast import Blast
from src.common.misc import check_evalue


class BlastN(Blast):

    def search_database(self, evalue=0.001):
        check_evalue(evalue)
        cline = NcbiblastnCommandline(query=self.query, db=self.database, evalue=evalue, outfmt=5, out=self.output)
        cline()

    def search_subject(self, evalue=0.001):
        check_evalue(evalue)
        cline = NcbiblastnCommandline(query=self.query, subject=self.subject, evalue=evalue, outfmt=5, out=self.output)
        cline()
