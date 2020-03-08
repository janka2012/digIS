from ..common.Classifier import Classifier
from ..hmmer.HmmerHspFlat import HmmerHspFlat

class RecordDigISAttrib:

    def __init__(self, source_type, hmmer_hsp):
        self.source_type = source_type
        self.hmmer_hsp = hmmer_hsp
        self.status = 'valid'
        self.extension_level = 'domain'
        self.classification = None

    @classmethod
    def get_csv_header(cls):
        attrib_header = ["source_type", "status", "extension_level"]
        hsp_header = HmmerHspFlat.get_csv_header()
        class_header = Classifier.get_csv_header(verbose=True)
        return attrib_header + hsp_header + class_header

    def to_csv(self):
        row = [self.source_type, self.status, self.extension_level]
        hsp_row = self.hmmer_hsp.to_csv()
        class_header = Classifier.get_csv_header(verbose=True)
        class_row = ['' for _ in range(len(class_header))]
        if self.classification:
            class_row = self.classification.to_csv(verbose=True)
        return row + hsp_row + class_row