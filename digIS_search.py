import argparse

from src.search_tool.digISConfiguration import digISConfiguration
from src.search_tool.digISMultifasta import digISMultifasta


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="digIS search")

    parser.add_argument('-i', "--input", action='store', dest='input_fasta', required=True,
                        help='Input fasta file, nucleotide')

    parser.add_argument('-g', "--genbank", action='store', dest='genbank_file', required=False, default=None,
                        help='Genbank annotations for genome in the input fasta file.')

    parser.add_argument('-o', "--output", action='store', dest='output_dir', required=False, default="digIS_output",
                        type=str, help='Output directory name, default=digIS_output.')

    parser.add_argument('-n', '--no-translate', dest='translate', required=False, action='store_false')
    parser.set_defaults(translate=True)

    parser.add_argument('-f', "--format", action='store', dest='out_format', required=False, default="csv", type=str,
                        choices=["csv", "gff"], help='Output format, default csv. Possible choices: csv, gff.')

    args = parser.parse_args()

    digIS_conf = digISConfiguration(genome_file=args.input_fasta,
                                    genbank_file=args.genbank_file,
                                    out_format=args.out_format,
                                    output_dir=args.output_dir)

    digIS = digISMultifasta(digIS_conf)
    digIS.run()

