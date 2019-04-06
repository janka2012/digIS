import argparse

from .src.search_tool.digIS import digIS
from .src.search_tool.digISConfiguration import digISConfiguration


def print_args():
    print('input fasta =', args.input_fasta)
    print('genbank file =', args.genbank_file)
    print('context size orf =', args.context_size_orf)
    print('context size is =', args.context_size_is)
    print('max merge distance =', args.merge_distance)
    print('min genbank overlap =', args.genbank_overlap)
    print('output dir =', args.output_dir)
    print('output fmt =', args.out_format)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="digIS search")

    parser.add_argument('-i', "--input", action='store', dest='input_fasta', required=True,
                        help='Input fasta file, nucleotide')

    parser.add_argument('-g', "--genbank", action='store', dest='genbank_file', required=False, default=None,
                        help='Genbank annotations for genome in the input fasta file.')

    parser.add_argument('-c', "--contextORF", action='store', dest='context_size_orf', required=False, default=1000,
                        type=int, help='Context size, default=1000 bp')

    parser.add_argument('-s', "--contextIS", action='store', dest='context_size_is', required=False, default=2000,
                        type=int, help='Context size, default=2000 bp')

    parser.add_argument('-d', "--merge_distance", action='store', dest='merge_distance', required=False, default=300,
                        type=int, help='Maximum distance for merging hits, default=300 bp')

    parser.add_argument('-b', "--genbank_overlap", action='store', dest='genbank_overlap', required=False, default=100,
                        type=int, help='Minimum overlap with genbank annotation, default=100 bp')

    parser.add_argument('-o', "--output", action='store', dest='output_dir', required=False, default="digIS_output",
                        type=str, help='Output directory name, default=digIS_output.')

    parser.add_argument('-f', "--format", action='store', dest='out_format', required=False, default="csv", type=str,
                        choices=["csv", "gff"], help='Output format, default csv. Possible choices: csv, gff.')

    args = parser.parse_args()

    digis_conf = digISConfiguration(genome_file=args.input_fasta,
                                    context_size_orf=args.context_size_orf,
                                    context_size_is=args.context_size_is,
                                    max_merge_distance=args.merge_distance,
                                    genbank_file=args.genbank_file,
                                    min_gb_overlap=args.genbank_overlap,
                                    out_format=args.out_format,
                                    output_dir=args.output_dir)

    dIS = digIS(digis_conf)
    dIS.run()
