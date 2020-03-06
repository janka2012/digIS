class HmmerHspFlat:
    
    def __init__(self, query, hit, hsp):
        self.sid = hit.id
        self.slen = hit.seq_len

        self.qid = query.id
        self.qlen = query.seq_len

        self.seq_evalue = float(hit.evalue)
        self.seq_bitscore = float(hit.bitscore)
        self.seq_bias = float(hit.bias)

        self.dom_idx = hsp.domain_index
        self.dom_num = len(hit.hsps)
        self.dom_cevalue = float(hsp.evalue_cond)
        self.dom_evalue = float(hsp.evalue)
        self.dom_bitscore = float(hsp.bitscore)
        self.dom_bias = float(hsp.bias)

        self.qstart = hsp.query_start + 1
        self.qend = hsp.query_end
        self.sstart = hsp.hit_start + 1
        self.send = hsp.hit_end
        self.sstart_env = hsp.env_start + 1
        self.send_env = hsp.env_end
        self.acc = float(hsp.acc_avg)
        self.sdesc = hit.description

    def to_csv(self):
        header = ["subject_id", "subject_len", "query_id", "query_len", "seq_evalue", "seq_bitscore", "seq_bias",
                  "dom_idx", "dom_num", "dom_cevalue", "dom_evalue", "dom_bitscore", "dom_bias", "query_start",
                  "query_end", "subject_start", "subject_end", "subject_env_start", "subject_env_end", "acc_avg",
                  "subject_description"]

        row = [self.sid, self.slen, self.qid, self.qlen, self.seq_evalue, self.seq_bitscore, self.seq_bias,
                self.dom_idx, self.dom_num, self.dom_cevalue, self.dom_evalue, self.dom_bitscore, self.dom_bias,
                self.qstart, self.qend, self.sstart, self.send, self.sstart_env, self.send_env, self.acc, self.sdesc]

        return header, row