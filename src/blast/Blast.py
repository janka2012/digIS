import os
import tempfile

from abc import ABC
from abc import abstractmethod
from copy import copy
from Bio import SeqIO
from Bio.Blast import NCBIXML

from src.blast.BlastHspFlat import BlastHspFlat
from src.blast.BlastQuery import BlastQuery


class Blast(ABC):

    def __init__(self, query, database="", subject="", output=None, remove_query=False):
        self.query = query
        self.subject = subject
        self.database = database
        self.best_hit = BlastHspFlat()
        self.best_hits = {}
        self.query_hits = []
        self.remove_query = remove_query

        if output:
            self.output = output
            self.remove_output = False
        else:
            self.fd, self.output = tempfile.mkstemp(prefix="digIS_", suffix=".xml")
            self.remove_output = True

    @abstractmethod
    def search_database(self):
        raise NotImplementedError("Should have implemented this")

    @abstractmethod
    def search_subject(self):
        raise NotImplementedError("Should have implemented this")

    @classmethod
    def from_seqrec(cls, qrec, database, output=None):
        fd, query = tempfile.mkstemp(prefix="digIS_", suffix=".fasta")
        SeqIO.write(qrec, query, "fasta")
        os.close(fd)
        return cls(query, database=database, output=output, remove_query=True)

    def parse(self):
        with open(self.output) as result_handle:
            try:
                blast_records = NCBIXML.parse(result_handle)
                self.query_hits = []
                for record in blast_records:
                    query_hits = BlastQuery.from_rec(record)
                    self.query_hits.append(query_hits)
            except ValueError:
                print("XML file with blast results is empty")

    def get_best_hit(self):
        for rec in self.query_hits:
            hspflat = rec.get_best_hit()
            if hspflat.score > 0.0:
                self.best_hits[rec.query_id] = copy(hspflat)
            if self.best_hit < hspflat:
                self.best_hit = hspflat
        return self.best_hit

    def print_best_hits(self):
        print('\n'.join(list(str(i) for i in self.best_hits.values())))

    def __del__(self):
        if self.remove_output and os.path.exists(self.output):
            os.close(self.fd)
            os.remove(self.output)
        if self.remove_query and os.path.exists(self.query):
            os.remove(self.query)

    def __str__(self):
        return '\n'.join(list(str(i) for i in self.query_hits))
