# digIS


**Table of content**
<!---toc start-->

  * [Requirements](#requirements)
  * [Installing](#installing)
  * [Usage](#usage)
  * [Understanding Outputs](#understanding-outputs)

<!---toc end-->

## Requirements
- HMMER 3.1b2 or higher, download the latest version from http://hmmer.org/download.html
- ncbi-blast+ v 2.6.1 or higher, download the latest version from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- Python3 (Biopython 1.73 or higher)

## Install dependencies using package manager (for Ubuntu)
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

## Download digIS from github repository
```bash
git clone https://github.com/janka2012/digIS.git
```

## Usage

### Mode with GenBank annotation

```bash
export PYTHONPATH=/path/to/digis/
python3 digIS_search.py -i data/test_data/NC_002608.fasta -g data/test_data/NC_002608.gb -o digis_genbank -f csv
```

### Mode without GenBank annotation
```bash
export PYTHONPATH=/path/to/digis/
python3 digIS_search.py -i data/test_data/NC_002608.fasta -o digis_genbank -f csv
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
sh digis_docker_wrapper.sh -i data/test_data/NC_002608.fasta -g data/test_data/NC_002608.gb -o digis_genbank -f csv
```

### Understanding Outputs

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
