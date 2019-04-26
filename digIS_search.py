import argparse

from src.search_tool.digIS import digIS
from src.search_tool.digISConfiguration import digISConfiguration


def print_args():
    print('input fasta =', args.input_fasta)
    print('genbank file =', args.genbank_file)
    print('min genbank overlap =', args.genbank_overlap)
    print('output dir =', args.output_dir)
    print('output fmt =', args.out_format)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="digIS search")

    parser.add_argument('-i', "--input", action='store', dest='input_fasta', required=True,
                        help='Input fasta file, nucleotide')

    parser.add_argument('-g', "--genbank", action='store', dest='genbank_file', required=False, default=None,
                        help='Genbank annotations for genome in the input fasta file.')

    parser.add_argument('-b', "--genbank_overlap", action='store', dest='genbank_overlap', required=False, default=100,
                        type=int, help='Minimum overlap with genbank annotation, default=100 bp')

    parser.add_argument('-o', "--output", action='store', dest='output_dir', required=False, default="digIS_output",
                        type=str, help='Output directory name, default=digIS_output.')

    parser.add_argument('-f', "--format", action='store', dest='out_format', required=False, default="csv", type=str,
                        choices=["csv", "gff"], help='Output format, default csv. Possible choices: csv, gff.')

    args = parser.parse_args()

    digis_conf = digISConfiguration(genome_file=args.input_fasta,
                                    genbank_file=args.genbank_file,
                                    min_gb_overlap=args.genbank_overlap,
                                    out_format=args.out_format,
                                    output_dir=args.output_dir)

    dIS = digIS(digis_conf)
    dIS.run()
