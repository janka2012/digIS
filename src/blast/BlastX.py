from Bio.Blast.Applications import NcbiblastxCommandline

from definitions import BLASTX_GAPOPEN, BLASTX_GAPEXTEND, BLASTX_WORDSIZE, BLASTX_EVALUE
from ..blast.Blast import Blast
from ..common.misc import check_evalue


class BlastX(Blast):

    def search_database(self, evalue=BLASTX_EVALUE):
        check_evalue(evalue)
        cline = NcbiblastxCommandline(query=self.query, db=self.database, evalue=evalue, outfmt=5, out=self.output,
                                      gapopen=BLASTX_GAPOPEN, gapextend=BLASTX_GAPEXTEND, word_size=BLASTX_WORDSIZE)
        cline()

    def search_subject(self, evalue=BLASTX_EVALUE):
        check_evalue(evalue)
        cline = NcbiblastxCommandline(query=self.query, subject=self.subject, evalue=evalue, outfmt=5, out=self.output,
                                      gapopen=BLASTX_GAPOPEN, gapextend=BLASTX_GAPEXTEND, word_size=BLASTX_WORDSIZE)
        cline()
