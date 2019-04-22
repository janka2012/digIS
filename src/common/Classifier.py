
class Classifier:

    def __init__(self, rec, gb_rec, orf_blast_hit, is_blast_hit):
        self.rec = rec
        self.genbank_recs = gb_rec
        self.blast_orf = orf_blast_hit
        self.blast_is_dna = is_blast_hit
        self.similarity_orf = None
        self.similarity_is = None
        self.similarity_all = None
        self.genbank_annotation = None
        self.level = None

    def classify(self):
        self.__assign_overall_similarity_with_isfinderdb()
        if self.genbank_recs:
            self.__clean_duplicit_gene_records()
            self.__assign_genbank_annotation()
            self.__assign_level()

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

    def __assign_genbank_annotation(self):
        if self.__no_genbank_annotation():
            genbank_annotation = 'no'
        elif self.__is_annotated_IS():
            genbank_annotation = 'is_related'
        elif self.__is_hypotetical_IS():
            genbank_annotation = 'no'
        else:
            genbank_annotation = 'other_record'
        self.genbank_annotation = genbank_annotation

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

    def __assign_overall_similarity_with_isfinderdb(self):
        self.similarity_is = self.__assign_similarity_level_dna()
        self.similarity_orf = self.__assign_similarity_level_orf()
        dict_orf_is_all = {('weak', 'weak'): 'weak', ('weak', 'medium'): 'medium', ('weak', 'strong'): 'strong',
                           ('medium', 'weak'): 'medium', ('medium', 'medium'): 'medium', ('medium', 'strong'): 'strong',
                           ('strong', 'weak'): 'strong', ('strong', 'medium'): 'strong', ('strong', 'strong'): 'strong'}

        self.similarity_all = dict_orf_is_all[(self.similarity_orf, self.similarity_is)]

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

    def __assign_level(self):
        dict_gb_sim_all = {('no', 'weak'): 'wFP',
                           ('no', 'medium'): 'pNov',
                           ('no', 'strong'): 'wTP',
                           ('is_related', 'weak'): 'wTP',
                           ('is_related', 'medium'): 'wTP',
                           ('is_related', 'strong'): 'sTP',
                           ('other_record', 'weak'): 'sFP',
                           ('other_record', 'medium'): 'wFP',
                           ('other_record', 'strong'): 'wTP'}

        self.level = dict_gb_sim_all[self.genbank_annotation, self.similarity_all]

    def to_csv(self):
        header = ["class_sim_orf", "class_sim_is", "class_sim_all", "class_genebank", "class_level"]
        row = [self.similarity_orf, self.similarity_is, self.similarity_all, self.genbank_annotation, self.level]
        return header, row