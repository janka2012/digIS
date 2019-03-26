class digISClassification:

    def __init__(self):
        self.similarity_orf = ""
        self.similarity_is = ""
        self.similarity_all = ""
        self.genbank_annotation = None
        self.level = tuple()

    def __str__(self):
        return "Sim ORF: {}, Sim IS: {}, Sim All: {}, GB annot: {}, Level: {}".format(self.similarity_orf,
                                                                                      self.similarity_is,
                                                                                      self.similarity_all,
                                                                                      self.genbank_annotation,
                                                                                      self.level)
