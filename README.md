# Gibberellin-Operon
## Purpose

ROAGUE is a tool to reconstruct ancestors of gene blocks in prokaryotic genomes. This repo is an application of [ROAGUE](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction) on operon Gibberellin. 

## Requirements
* [Wget](https://www.gnu.org/software/wget/) 
* [Conda](https://conda.io/miniconda.html) (package manager so we don't have to use sudo)
* [Python 3+](https://www.python.org/download/releases/3.0/)
* [Biopython 1.63+](http://biopython.org/wiki/Download)
* [Clustalw](http://www.clustal.org/clustal2/#Download)
* [Muscle Alignment](https://www.drive5.com/muscle/downloads.htm)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [ETE3](http://etetoolkit.org/download/) (python framework for tree)
* [PDA](http://www.cibiv.at/software/pda/#download) (optional if you want to debias your tree base on Phylogenetic Diversity)
* [Prokka](http://www.vicbioinformatics.com/software.shtml) (bacterial genome annotation tool)

## Installation
Users can either use github interface Download button or type the following command in command line:
```bash
git clone https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction
```
Install Miniconda (you can either export the path everytime you use ROAGUE, or add it to the .bashrc file). Before using
the following command line, users will need to install [Wget](https://www.gnu.org/software/wget/).
```bash
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O Miniconda-latest-Linux-x86_64.sh
bash Miniconda-latest-Linux-x86_64.sh -b -p ~/anaconda_ete/
export PATH=~/anaconda_ete/bin:$PATH;

```

Install Biopython and ete3 using conda (highly recommended install biopython with conda)
```bash
conda install -c bioconda biopython ete3
```
Install ete_toolchain for visualization
```bash
conda install -c etetoolkit ete_toolchain
```

Install BLAST, ClustalW, MUSCLE, Prokka
```bash
conda install -c bioconda blast clustalw muscle prokka
```

For PDA, check installation instructions on this website: [PDA](http://www.cibiv.at/software/pda/#download)

## Usage

The easiest way to run the project is to execute the script [ROAGUE](https://github.com/nguyenngochuy91/Gibberellin-Operon/blob/master/ROAGUE.py), 
the results are stored in directory **ancestral_reconstruction**, under the name **Ryan_operon**, **Ryan_rpoB**,**Ryan_operon_64**,**Ryan_rpoB_64**
```ba
./ROAGUE.py
```
