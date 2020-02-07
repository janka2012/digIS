from Bio.Blast.Applications import NcbiblastnCommandline

from definitions import BLASTN_GAPOPEN, BLASTN_GAPEXTEND, BLASTN_WORDSIZE, BLASTN_EVALUE
from ..blast.Blast import Blast
from ..common.misc import check_evalue


class BlastN(Blast):

    def search_database(self, evalue=BLASTN_EVALUE):
        check_evalue(evalue)
        cline = NcbiblastnCommandline(query=self.query, db=self.database, evalue=evalue, outfmt=5, out=self.output,
                                      gapopen=BLASTN_GAPOPEN, gapextend=BLASTN_GAPEXTEND, word_size=BLASTN_WORDSIZE)
        cline()

    def search_subject(self, evalue=BLASTN_EVALUE):
        check_evalue(evalue)
        cline = NcbiblastnCommandline(query=self.query, subject=self.subject, evalue=evalue, outfmt=5, out=self.output,
                                      gapopen=BLASTN_GAPOPEN, gapextend=BLASTN_GAPEXTEND, word_size=BLASTN_WORDSIZE)
        cline()
