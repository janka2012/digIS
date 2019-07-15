import argparse

from collections import OrderedDict

from src.search_tool.digIS import digIS
from src.search_tool.digISConfiguration import digISConfiguration
from src.common.genbank import read_gb
from src.common.genome import Genome


def print_args():
    print('input fasta =', args.input_fasta)
    print('genbank file =', args.genbank_file)
    print('output dir =', args.output_dir)
    print('output fmt =', args.out_format)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="digIS search")

    parser.add_argument('-i', "--input", action='store', dest='input_fasta', required=True,
                        help='Input fasta file, nucleotide')

    parser.add_argument('-g', "--genbank", action='store', dest='genbank_file', required=False, default=None,
                        help='Genbank annotations for genome in the input fasta file.')

    parser.add_argument('-o', "--output", action='store', dest='output_dir', required=False, default="digIS_output",
                        type=str, help='Output directory name, default=digIS_output.')

    parser.add_argument('-t', '--translate', dest='translate', required=False, action='store_true')
    parser.set_defaults(translate=False)

    parser.add_argument('-c', '--classify', dest='classify', required=False, action='store_true')
    parser.set_defaults(classify=False)

    parser.add_argument('-f', "--format", action='store', dest='out_format', required=False, default="csv", type=str,
                        choices=["csv", "gff"], help='Output format, default csv. Possible choices: csv, gff.')

    args = parser.parse_args()

    digIS_conf = digISConfiguration(genome_file=args.input_fasta,
                                    genbank_file=args.genbank_file,
                                    out_format=args.out_format,
                                    output_dir=args.output_dir,
                                    translate=args.translate)

    genomes_dict = Genome.parse_genomes(fasta_file=digIS_conf.genome_file)
    genbank_dict = read_gb(digIS_conf.genbank_file) if digIS_conf.genbank_file else OrderedDict()

    for i, (genome_id, genome_rec) in enumerate(genomes_dict.items(), start=1):
        dIS = digIS(digIS_conf, genome=genome_rec, genbank_features=genbank_dict.get(genome_id, []))
        dIS.run(search=True, classify=True)

