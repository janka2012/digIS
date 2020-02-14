from Bio.Blast.Applications import NcbiblastxCommandline

from definitions import BLASTX_GAPOPEN, BLASTX_GAPEXTEND, BLASTX_WORDSIZE, BLASTX_EVALUE, NUM_THREADS
from ..blast.Blast import Blast
from ..common.misc import check_evalue


class BlastX(Blast):

    def __get_default_params(self):
        params = {'query': self.query, 'outfmt': 5, 'out': self.output, 'gapopen': BLASTX_GAPOPEN,
                  'gapextend': BLASTX_GAPEXTEND, 'word_size': BLASTX_WORDSIZE}

        if NUM_THREADS > 0:
            params['num_threads'] = NUM_THREADS

        return params

    def search_database(self, evalue=BLASTX_EVALUE):
        check_evalue(evalue)

        params = self.__get_default_params()
        params['db'] = self.database
        params['evalue'] = evalue

        cline = NcbiblastxCommandline(**params)
        cline()

    def search_subject(self, evalue=BLASTX_EVALUE):
        check_evalue(evalue)

        params = self.__get_default_params()
        params['subject'] = self.subject
        params['evalue'] = evalue

        cline = NcbiblastxCommandline(**params)
        cline()
