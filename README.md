# digIS

## Requirements
- HMMER 3.1b2 or higher, download the latest version from http://hmmer.org/download.html
- ncbi-blast+ v 2.8.1 or higher, download the latest version from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- Python3 (Biopython 1.73 or higher)

## Installing
```
sudo apt-get update
sudo apt-get install hmmer
sudo apt-get install ncbi-blast+
```

## Searching

### Mode with GenBank annotation

### Mode without GenBank annotation

## Understanding Ouputs

### Getting FASTA file using GFF file

The [GFF](http://gmod.org/wiki/GFF3) is a standard format for storing of genome features. This file can be used as an input for other tools to process or visualize appropriate genomic features. 

For instance, FASTA sequences of IS elements (their catalytic domains) can be obtained using [Bedtools](https://bedtools.readthedocs.io/en/latest/) and command `getfasta` as follows:     

```
bedtools getfasta -fi <input.fasta> -bed <input.gff> -fo <output.fasta>
```
where _input.fasta_ represents FASTA file used for searching, _input.gff_ is the digIS output GFF file and _output.fasta_ is required output file. 

### Getting Flank Regions using GFF file

```
bedtools flank -i <input.gff> -g <genome.size> -l <left flank size>  -r <right flank size> 
```
where _genome.size_ is a text file containing information about chromosomes and their sizes in form: `chromosome_name<TAB>chromosome_size`. For more information about _genome.size_ file format please see Bedtools [documentation](https://bedtools.readthedocs.io/en/latest/).

As the output a new GFF file with positions of flank regions is generated. Then, the appropriate FASTA file with flank sequences can be obtained using `bedtools getfasta` command described above. 
