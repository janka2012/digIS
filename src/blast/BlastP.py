from Bio.Blast.Applications import NcbiblastpCommandline
from ..blast.Blast import Blast
from ..common.misc import check_evalue


class BlastP(Blast):

    def search_database(self, evalue=0.001):
        check_evalue(evalue)
        try:
            cline = NcbiblastpCommandline(query=self.query, db=self.database, evalue=evalue, outfmt=5, out=self.output)
            cline()
        except ValueError:
            print("query file is empty.")

    def search_subject(self, evalue=0.001):
        check_evalue(evalue)
        cline = NcbiblastpCommandline(query=self.query, subject=self.subject, evalue=evalue, outfmt=5, out=self.output)
        cline()
