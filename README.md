# digIS


**Table of content**
<!---toc start-->

  * [Overview](#overview)
  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Understanding Outputs](#understanding-outputs)

<!---toc end-->

## Overview
digIS is a command-line tool for detection of insertion sequences in prokaryotic genomes. It was developed in Python3 and utilizes several external tools such as HMMER v3.1, BLAST, and Biopython library. As an input, digIS accepts genomic sequences (e.g. contigs, fully assembled prokaryotic genomes or other DNA sequences) in the FASTA format. Optionally, the user can also provide a GenBank annotation file for a given input sequence(s). This annotation is later used to improve classification of identified IS elements.

digIS search pipeline operates in following steps: (1) the whole input DNA sequence is translated into protein sequences (all possible six frames), (2) the protein sequences are searched using manually curated pHMM models and thus the seeds of putative IS elements are identified, (3) the seeds are extended according to sequence similarity with known IS elements in ISFInder database, (4) extended seeds are classified based on sequence similarity and optionally using GenBank annotations to help assess their quality, (5) finally, the classified outputs are reported in CSV and GFF3 format.

## Requirements
- HMMER 3.1b2 or higher, download the latest version from http://hmmer.org/download.html
- ncbi-blast+ v 2.6.1 or higher, download the latest version from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- Python3 (Biopython 1.73 or higher)

## Installation

### Install dependencies using package manager (for Ubuntu)
```bash
sudo apt-get update
sudo apt-get install hmmer
sudo apt-get install ncbi-blast+

# install python3  and pip3
sudo apt-get install python3.7
sudo apt install python3-pip

# download Biopython
pip3 install biopython
```

### Download digIS version from github repository
```bash
# download the latest version
git clone https://github.com/janka2012/digIS.git

# or download specific release
wget https://github.com/janka2012/digIS/archive/v1.0.tar.gz
tar -xvzf v1.0.tar.gz
```
## Usage

### Mode with GenBank annotation

```bash
export PYTHONPATH=$PYTHONPATH:/path/to/digis/
python3 digIS_search.py -i data/test_data/NC_002608.fasta -g data/test_data/NC_002608.gb -o digis_genbank
```

### Mode without GenBank annotation
```bash
export PYTHONPATH=$PYTHONPATH:/path/to/digis/
python3 digIS_search.py -i data/test_data/NC_002608.fasta -o digis_without_genbank
```

## Run digIS in docker container

### Install docker

```bash
# update software repositories
sudo apt-get update

# uninstall older versions of docker
sudo apt-get remove docker docker-engine docker.io

# install docker
sudo apt install docker.io

# start and automate docker
sudo systemctl start docker
sudo systemctl enable docker

# check docker version (optional)
docker --version
```

### Build a digIS docker image

```bash
cd /path/to/digIS
docker build -t digis .

# List created docker images. You should see image with name digis listed.
docker images
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
digis               latest              15ef3194ff7d        27 hours ago        764MB
```

### Run digIS using ```digis_docker_wrapper.sh```
Instead of typing overwhelmingly long docker commands we are providing `digis_docker_wrapper.sh` script which allows you to use digIS docker image in really convinient way. This script takes same arguments as standard `digIS.py` script.

```bash
sh digis_docker_wrapper.sh -i data/test_data/NC_002608.fasta -g data/test_data/NC_002608.gb -o digis_genbank
```

## Understanding Outputs

### digIS output directory structure
digIS stores all results in output directory you specify by using `-o` option. The default output directory name is set to `digIS_output`. The output directory has following structure:
* `pep`: translated protein sequences of the input genomic sequences
* `hmmer`: results from hmmer and phmmer search, separate file for each sequence in input (multi)fasta file 
* `logs`: 
* `results`: digIS results for all sequences in input (multi)fasta file in CSV and GFF3 format together with summarization statistics.
    * `.csv`:
    * `.gff`:
    * `.sum`:

### digIS output files in `results` subfolder
For a given input (multi)fasta file digIS generates three files with results: CSV file, GFF3 file and file with summary statistics.

#### CSV output

* `id`:
* `level`: 
* `qid`: name of the query profile
* `qstart`: start position of the hit in the query profile
* `qend`: end position of the hit in the query profile
* `sid`: name of the subject sequence
* `sstart`: start position of the hit in subject sequence
* `send`: end position of the hit in subject sequence
* `strand`: the strand on which the hit was found
* `acc`: taken from Hmmer. A measure of how reliable the overall alignment is (from 0 to 1, with 1.00 indicating a completely reliable alignment according to the model).
* `class_sim_orf`: classification based on protein sequence similarity. Possible values: strong, medium, weak.
* `class_sim_is`: classification based on nucleotide sequence similarity. Possible values: strong, medium, weak.
* `class_sim_all`: overall classification, maximum from class_sim_orf and class_sim_is. Possible values: strong, medium, weak.
* `class_genebank`: classification based on GenBank annotation. Possible values: is_related, no, other record. If GenBank annotation is not provided, this classification is not available and this field is empty.
* `class_level`: overall classification resulting from class_sim_all and class_genebank. Possible values: sTP, wTP, pNov, wFP, wFP. If GenBank annotation is not provided, this classification is not available and this field is empty.

#### GFF3 output
* seqid
* source
* type
* start
* end
* score
* strand
* frame
* attribute

#### Summary statistics


## Getting FASTA file using GFF file

The [GFF](http://gmod.org/wiki/GFF3) is a standard format for storing of genome features. This file can be used as an input for other tools to process or visualize appropriate genomic features. 

For instance, FASTA sequences of IS elements (their catalytic domains) can be obtained using [Bedtools](https://bedtools.readthedocs.io/en/latest/) and command `getfasta` as follows:     

```bash
bedtools getfasta -fi <input.fasta> -bed <input.gff> -fo <output.fasta>
```
where _input.fasta_ represents FASTA file used for searching, _input.gff_ is the digIS output GFF file and _output.fasta_ is required output file. 

### Getting Flank Regions using GFF file

```bash
bedtools flank -i <input.gff> -g <genome.size> -l <left flank size>  -r <right flank size> 
```
where _genome.size_ is a text file containing information about chromosomes and their sizes in form: `chromosome_name<TAB>chromosome_size`. For more information about _genome.size_ file format please see Bedtools [documentation](https://bedtools.readthedocs.io/en/latest/).

As the output a new GFF file with positions of flank regions is generated. Then, the appropriate FASTA file with flank sequences can be obtained using `bedtools getfasta` command described above. 
