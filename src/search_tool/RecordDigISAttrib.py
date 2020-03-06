from ..common.Classifier import Classifier

class RecordDigISAttrib:

    def __init__(self, source_type, hmmer_hsp):
        self.source_type = source_type
        self.hmmer_hsp = hmmer_hsp
        self.status = 'valid'
        self.extension_level = 'domain'
        self.classification = None

    def to_csv(self):
        header = ["source_type", "status", "extension_level"]
        row = [self.source_type, self.status, self.extension_level]
        hsp_header, hsp_row = self.hmmer_hsp.to_csv()
        class_header = Classifier.get_csv_header(verbose=True)
        class_row = ['' for _ in range(len(class_header))]
        if self.classification:
            class_header, class_row = self.classification.to_csv(verbose=True)
        return header + hsp_header + class_header, row + hsp_row + class_row