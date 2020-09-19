import re

from definitions import IS_GB_KEYWORDS, HYPOTHETICAL_GB_KEYWORDS, IS_FAMILIES_NAMES, INTRAFAMILY_ORF_SIM_THRESHOLD, INTERFAMILY_ORF_SIM_THRESHOLD


class Classifier:

    def __init__(self, rec, gb_rec, orf_blast_hit, is_blast_hit):
        self.rec = rec
        self.genbank_recs = gb_rec
        self.blast_orf = orf_blast_hit
        self.blast_is_dna = is_blast_hit
        self.similarity_orf = None
        self.similarity_is = None
        self.genbank_annotation = None
        self.level = None
        self.kept = True

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

            if "incomplete" in ",".join(gb_annots):
                out = False
            elif (rec.type in ['mobile_element', 'mobile_element_type'] or rec.type in ['CDS', 'gene', 'misc_feature']) \
                    and any(annot.lower() in ",".join(gb_annots) for annot in IS_GB_KEYWORDS + IS_FAMILIES_NAMES):
                out = True
            elif (rec.type in ['mobile_element', 'mobile_element_type'] or rec.type in ['CDS', 'gene', 'misc_feature']) \
                    and 'integrase' in ",".join(gb_annots) and self.blast_orf.subject_identity >= 0.9:
                out = True
        return out

    def __is_hypotetical_IS(self):
        hp_all_length = 0
        other_all_length = 0
        for rec in self.genbank_recs:
            gb_annots_product, gb_annots_note = self.__get_genbank_annotations(rec.qualifiers)
            gb_annots = gb_annots_product + gb_annots_note
            overlap_length = rec.get_overlap_length(self.rec)
            if rec.type in ['CDS', 'gene'] and any(annot.lower() in ",".join(gb_annots) for annot in HYPOTHETICAL_GB_KEYWORDS):
                hp_all_length += overlap_length
            elif re.match(r'DUF\d+', ",".join(gb_annots_product), re.M | re.I):
                hp_all_length += overlap_length
            elif rec.type == 'CDS' and not gb_annots_product and gb_annots_note:
                hp_all_length += overlap_length
            elif "incomplete" in ",".join(gb_annots):
                other_all_length += overlap_length
            else:
                other_all_length += overlap_length

        hp_all_coverage = 0
        if hp_all_length > 0:
            hp_all_coverage = self.rec.width/hp_all_length

        out = True if hp_all_coverage >= 0.50 else False
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

    def __assign_similarity_level_dna(self):
        return self.blast_is_dna.subject_identity

    def __assign_similarity_level_orf(self):
        return self.blast_orf.subject_identity

    def __assign_level(self):
        if self.genbank_annotation == "is_related":
            self.level = "TP"
        elif self.genbank_annotation == "other_record":
            self.level = "FP"
        else:
            if self.similarity_orf > INTRAFAMILY_ORF_SIM_THRESHOLD:
                self.level = "TP"
            elif self.similarity_orf < INTERFAMILY_ORF_SIM_THRESHOLD:
                self.level = "FP"
            else:
                self.level = "pNov"

    @classmethod
    def get_csv_header(cls, verbose=False):
        if verbose:
            header = ['Genome', 'Annotation', 'Orf_Sim', 'IS_Sim', 'Str_Rec', 'Str_GB', 'Str_Orf', 'Str_IS', 'kept']
        else:
            header = ["ORF_sim", "IS_sim", "GenBank_class"]
        return header

    def to_csv(self, verbose=False):
        if verbose:
            str_gb = ""
            if self.genbank_recs is not None:
                str_gb = '[' + ','.join(str(i) for i in self.genbank_recs) + ']' if len(self.genbank_recs) > 0 else ""
            str_bl_orf = str(self.blast_orf) if self.blast_orf.score != 0.0 else ""
            str_bl_is = str(self.blast_is_dna) if self.blast_is_dna.score != 0.0 else ""
            row = [self.rec.genome_name, self.genbank_annotation, self.similarity_orf,
                   self.similarity_is, str(self.rec), str_gb, str_bl_orf, str_bl_is, self.kept]
        else:
            row = [self.similarity_orf, self.similarity_is, self.genbank_annotation]
        return row

    def __eq__(self, other):
        return self.level == other.level

    def __lt__(self, other):
        return ['FP', 'pNov', 'TP'].index(self.level) < ['FP', 'pNov', 'TP'].index(other.level)

    def __gt__(self, other):
        return ['FP', 'pNov',  'TP'].index(self.level) > ['FP', 'pNov', 'TP'].index(other.level)

