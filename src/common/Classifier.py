import re

from definitions import IS_GB_KEYWORDS, HYPOTHETICAL_GB_KEYWORDS, IS_FAMILIES_NAMES


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
        if self.genbank_recs is not None:
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
        out = True
        if len(self.genbank_recs) == 0:
            return out

        for rec in self.genbank_recs:
            gb_annots_product, gb_annots_note = self.__get_genbank_annotations(rec.qualifiers)
            gb_annots = gb_annots_product + gb_annots_note

            if "functional annotations will be submitted" not in ",".join(gb_annots):
                out = False
                break
        return out

    def __is_annotated_IS(self):
        out = False

        for rec in self.genbank_recs:
            gb_annots_product, gb_annots_note = self.__get_genbank_annotations(rec.qualifiers)
            gb_annots = gb_annots_product + gb_annots_note

            if rec.type in ['mobile_element', 'mobile_element_type'] or \
                    rec.type in ['CDS', 'gene', 'misc_feature'] and \
                    any(annot.lower() in ",".join(gb_annots) for annot in IS_GB_KEYWORDS + IS_FAMILIES_NAMES):
                out = True
        return out

    def __is_hypotetical_IS(self):
        out = False

        for rec in self.genbank_recs:
            gb_annots_product, gb_annots_note = self.__get_genbank_annotations(rec.qualifiers)
            gb_annots = gb_annots_product + gb_annots_note

            if rec.type in ['CDS', 'gene'] and any(annot.lower() in ",".join(gb_annots) for annot in HYPOTHETICAL_GB_KEYWORDS):
                out = True
            elif re.match(r'DUF\d+', ",".join(gb_annots_product), re.M | re.I):
                out = True
            elif rec.type == 'CDS' and not gb_annots_product and gb_annots_note:
                out = True
            else:
                out = False
                break
        return out

    def __get_genbank_annotations(self, gb_qualifiers):
        gb_annots_product = []
        gb_annots_note = []
        if 'product' in gb_qualifiers:
            gb_annots_product = list(map(str.lower, gb_qualifiers['product']))
        if 'note' in gb_qualifiers:
            gb_annots_note = list(map(str.lower, gb_qualifiers['note']))

        return gb_annots_product, gb_annots_note

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
        elif self.blast_is_dna.subject_identity < 0.7:
            similarity_is = 'medium'
        else:
            similarity_is = 'strong'
        return similarity_is

    def __assign_similarity_level_orf(self):
        if self.blast_orf.subject_identity < 0.25:
            similarity_orf = 'weak'
        elif self.blast_orf.subject_identity < 0.45:
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

    @classmethod
    def get_csv_header(cls, verbose=False):
        if verbose:
            header = ['Genome', 'Level', 'Similarity', 'Annotation', 'Orf_Sim', 'IS_Sim', 'Str_Rec', 'Str_GB', 'Str_Orf', 'Str_IS']
        else:
            header = ["class_sim_orf", "class_sim_is", "class_sim_all", "class_genebank", "class_level"]
        return header

    def to_csv(self, verbose=False):
        if verbose:
            str_gb = ""
            if self.genbank_recs is not None:
                str_gb = '[' + ','.join(str(i) for i in self.genbank_recs) + ']' if len(self.genbank_recs) > 0 else ""
            str_bl_orf = str(self.blast_orf) if self.blast_orf.score != 0.0 else ""
            str_bl_is = str(self.blast_is_dna) if self.blast_is_dna.score != 0.0 else ""
            row = [self.rec.genome_name, self.level, self.similarity_all, self.genbank_annotation, self.similarity_orf,
                   self.similarity_is, str(self.rec), str_gb, str_bl_orf, str_bl_is]
        else:
            row = [self.similarity_orf, self.similarity_is, self.similarity_all, self.genbank_annotation, self.level]
        return row

    def __eq__(self, other):
        return self.level == other.level

    def __lt__(self, other):
        return ['sFP','wFP','pNov','wTP','sTP'].index(self.level) < ['sFP','wFP','pNov','wTP','sTP'].index(other.level)
