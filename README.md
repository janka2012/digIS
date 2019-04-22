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
- ncbi-blast+ v 2.8.1 or higher, download the latest version from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- Python3 (Biopython 1.73 or higher)

## Install dependencies using package manager (for Ubuntu)
```
sudo apt-get update
sudo apt-get install hmmer
sudo apt-get install ncbi-blast+

# install python3  and pip3
sudo apt-get install python3.7
sudo apt install python3-pip

# download Biopython
pip3 install biopython
```

## Install dependencies from source

Create a directory where you will download and install required software
```
mkdir -p $HOME/bin
cd $HOME/bin
```

### Download and compile HMMER (v3.1b2) from source

```
# download
wget http://eddylab.org/software/hmmer/hmmer-3.1b2.tar.gz
tar zxf hmmer-3.1b2.tar.gz
cd hmmer-3.1b2

# compile
./configure --prefix $HOME/bin/hmmer-3.1b2
make
make check  # optional step
make install
```

Add HMMER to your ```$PATH``` and ```~/.profile```

```
HMMER=$HOME/bin/hmmer-3.1b2/bin
echo $"export PATH=\$PATH:$HMMER" >> ~/.profile
```

### Download an compile ncbi-blast+ from source

```
cd $HOME/bin

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz
tar zxf ncbi-blast-2.9.0+-x64-linux.tar.gz

BLAST=$HOME/bin/ncbi-blast-2.9.0+/bin
echo $"export PATH=\$PATH:$BLAST" >> ~/.profile
```

Installing Python3

```
wget https://www.python.org/ftp/python/3.7.3/Python-3.7.3.tgz
tar xf Python-3.7.3.tgz

cd Python-3.7.3.tgz


./configure --prefix $HOME/bin/Python-3.7.3
make
make test # optional

# problem with zlib
make install
```

Building Python

```
https://www.systutorials.com/topic/how-to-install-the-zlib-library-in-ubuntu/
```


## Usage

### Mode with GenBank annotation

### Mode without GenBank annotation

## Understanding Outputs

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
