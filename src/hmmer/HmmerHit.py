from src.hmmer.HmmerHsp import HmmerHsp
from src.common.sequence import transform_range


class HmmerHit:
    def __init__(self, hit, query_len):
        self.query_id = hit.query_id
        self.subject_id = hit.id
        self.subject_desc = hit.description
        self.query_len = query_len
        self.subject_len = hit.seq_len
        self.bitscore = float(hit.bitscore)
        self.evalue = hit.evalue
        self.hsps = [HmmerHsp(hsp, self.query_len) for hsp in hit.hsps]

    def get_best_hsp_in_range(self, seqlen, query_range=(0, 0), min_overlap=1):

        best_hsp = None

        for hsp in self.hsps:

            if query_range == (0, 0):
                if best_hsp < hsp:
                    best_hsp = hsp
            else:
                hit_range = (hsp.qstart, hsp.qend)
                frame = int(hsp.qid.strip()[-1])
                hsp_dna_range = transform_range(hit_range[0], hit_range[1], frame, seqlen)
                max_start = max(hsp_dna_range[0], query_range[0])
                min_end = min(hsp_dna_range[1], query_range[1])
                overlap_len = min_end - max_start + 1

                # init
                if not best_hsp and overlap_len >= min_overlap:
                    best_hsp = hsp
                else:
                    if best_hsp:
                        if best_hsp < hsp and overlap_len >= min_overlap:
                            best_hsp = hsp

        return best_hsp

    def __str__(self):
        return '{}, {}, {}, {}, {}, {}'.format(self.query_id,
                                               self.subject_id,
                                               self.subject_len,
                                               self.bitscore,
                                               self.evalue,
                                               self.subject_desc)
