class digISClassifier:

    def __init__(self):
        self.genbank_recs = []
        self.blast_orf = None
        self.blast_is_dna = None

    def classify(self, digis_rec, gb_rec, orf_blast_hit, is_blast_hit):

        self.genbank_recs = gb_rec
        self.blast_orf = orf_blast_hit
        self.blast_is_dna = is_blast_hit
        self.__assign_overall_similarity_with_isfinderdb(digis_rec)

        if self.genbank_recs:
            self.__clean_duplicit_gene_records()
            self.__assign_genbank_annotation(digis_rec)
            self.__assign_level(digis_rec)

    def __clean_duplicit_gene_records(self):
        """
        Clean records that have gene and CDS at the same positions
        :return:
        """

        out_idx = []
        for i, ref_rec in enumerate(self.genbank_recs):
            discard = False
            if ref_rec.type == 'gene':
                for rec in self.genbank_recs:
                    if rec.type == 'CDS' and rec.start == ref_rec.start and rec.end == ref_rec.end:
                        discard = True
            if not discard:
                out_idx.append(i)

        new_recs = list(self.genbank_recs[i] for i in out_idx)
        self.genbank_recs = new_recs

    def __assign_genbank_annotation(self, digis_rec):
        if self.__no_genbank_annotation():
            genbank_annotation = 'no'
        elif self.__is_annotated_IS():
            genbank_annotation = 'is related'
        elif self.__is_hypotetical_IS():
            genbank_annotation = 'no'
        else:
            genbank_annotation = 'other record'
        digis_rec.classification.genbank_annotation = genbank_annotation

    def __no_genbank_annotation(self):
        if len(self.genbank_recs) == 0:
            return True
        return False

    def __is_annotated_IS(self):
        out = False
        for rec in self.genbank_recs:
            if rec.type in ['mobile_element'] \
                    or (rec.type == 'CDS' and 'transposase' in ",".join(rec.qualifiers['product'])) \
                    or (rec.type == 'CDS' and 'insertion element protein' in ",".join(rec.qualifiers['product'])) \
                    or (rec.type == 'CDS' and 'serine recombinase' in ",".join(rec.qualifiers['product'])):
                out = True
        return out

    def __is_hypotetical_IS(self):
        out = False
        for rec in self.genbank_recs:
            if rec.type == 'CDS' and 'hypothetical protein' in ",".join(rec.qualifiers['product']):
                out = True
            else:
                out = False
                break
        return out

    def __assign_overall_similarity_with_isfinderdb(self, digis_rec):
        digis_rec.classification.similarity_is = self.__assign_similarity_level_dna()
        digis_rec.classification.similarity_orf = self.__assign_similarity_level_orf()
        dict_orf_is_all = {('weak', 'weak'): 'weak', ('weak', 'medium'): 'medium', ('weak', 'strong'): 'strong',
                           ('medium', 'weak'): 'medium', ('medium', 'medium'): 'medium', ('medium', 'strong'): 'strong',
                           ('strong', 'weak'): 'strong', ('strong', 'medium'): 'strong', ('strong', 'strong'): 'strong'}

        digis_rec.classification.similarity_all = dict_orf_is_all[(digis_rec.classification.similarity_orf,
                                                                   digis_rec.classification.similarity_is)]

    def __assign_similarity_level_dna(self):
        if self.blast_is_dna.subject_identity < 0.5:
            similarity_is = 'weak'
        elif self.blast_is_dna.subject_identity < 0.8:
            similarity_is = 'medium'
        else:
            similarity_is = 'strong'
        return similarity_is

    def __assign_similarity_level_orf(self):
        if self.blast_orf.subject_identity < 0.3:
            similarity_orf = 'weak'
        elif self.blast_orf.subject_identity < 0.6:
            similarity_orf = 'medium'
        else:
            similarity_orf = 'strong'
        return similarity_orf

    def __assign_level(self, digis_rec):
        dict_gb_sim_all = {('no', 'weak'): 'wFP',
                           ('no', 'medium'): 'pNov',
                           ('no', 'strong'): 'wTP',
                           ('is related', 'weak'): 'wTP',
                           ('is related', 'medium'): 'wTP',
                           ('is related', 'strong'): 'sTP',
                           ('other record', 'weak'): 'sFP',
                           ('other record', 'medium'): 'wFP',
                           ('other record', 'strong'): 'wTP'}

        digis_rec.classification.level = dict_gb_sim_all[digis_rec.classification.genbank_annotation,
                                                         digis_rec.classification.similarity_all]
