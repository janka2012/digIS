from Bio.Blast.Applications import NcbiblastnCommandline

from definitions import BLASTN_GAPOPEN, BLASTN_GAPEXTEND, BLASTN_WORDSIZE, BLASTN_EVALUE, NUM_THREADS
from ..blast.Blast import Blast
from ..common.misc import check_evalue


class BlastN(Blast):

    def __get_default_params(self):
        params = {'query' : self.query, 'outfmt' : 5, 'out' : self.output, 'gapopen' : BLASTN_GAPOPEN,
                  'gapextend' : BLASTN_GAPEXTEND, 'word_size' : BLASTN_WORDSIZE}

        if NUM_THREADS > 0:
            params['num_threads'] = NUM_THREADS

        return params

    def search_database(self, evalue=BLASTN_EVALUE):
        check_evalue(evalue)

        params = self.__get_default_params()
        params['db'] = self.database
        params['evalue'] = evalue

        cline = NcbiblastnCommandline(**params)
        cline()

    def search_subject(self, evalue=BLASTN_EVALUE):
        check_evalue(evalue)

        params = self.__get_default_params()
        params['subject'] = self.subject
        params['evalue'] = evalue

        cline = NcbiblastnCommandline(**params)
        cline()
