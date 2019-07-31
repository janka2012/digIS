import os
from collections import OrderedDict

from ..common.csv_utils import write_csv
from ..common.genbank import read_gb
from ..common.genome import Genome
from ..common.gff_utils import write_gff
from ..search_tool.digIS import digIS


class digISMultifasta:

    def __init__(self, config):
        self.config = config
        self.genomes_dict = Genome.parse_genomes(fasta_file=config.genome_file, output_dir=config.output_dir)
        self.genbank_dict = read_gb(config.genbank_file) if config.genbank_file else OrderedDict()
        self.digIS_recs = OrderedDict()
        for genome_id, genome_rec in self.genomes_dict.items():
            self.digIS_recs[genome_id] = digIS(self.config, genome=genome_rec, genbank_features=self.genbank_dict.get(genome_id, []))

    def run(self):
        for genome_id, digIS_rec in self.digIS_recs.items():
            digIS_rec.run()
        records = self.export()
        return records

    def export(self):
        print("===== Exporting outputs =====")
        fasta_basename = os.path.splitext(os.path.basename(self.config.genome_file))[0]
        csv_header = []
        csv_rows = []
        for genome_id, digIS_rec in self.digIS_recs.items():
            csv_header, rows = digIS_rec.export_records()
            csv_rows.extend(rows)

        print("Exporting records...")
        output_recs = os.path.join(self.config.output_dir, "results", fasta_basename + "." + self.config.out_format)
        if self.config.out_format == "csv":
            write_csv(csv_rows, output_recs, csv_header)
        elif self.config.out_format == "gff":
            write_gff(csv_rows, output_recs, csv_header)

        print("Exporting summary statistics...")
        sum_recs = []
        for genome_id, digIS_rec in self.digIS_recs.items():
            sum_recs.extend(digIS_rec.export_summary_stats())

        output_sum = os.path.join(self.config.output_dir, "results", fasta_basename + ".sum")
        with open(output_sum, 'w+', newline='') as f:
            f.write("{:40} {:15} {:>5} {:>10} {:>10} {:>10}\n".format('#seqid', 'family', 'nIS', 'bps', 'dnaLen', '%dna'))
            for rec in sum_recs:
                f.write("{:40} {:15} {:>5} {:>10} {:>10} {:>10.2f}\n".format(*rec))

        return len(csv_rows)
