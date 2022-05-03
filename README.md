# deepStrain: #

**deepStrain** is a strain-level metagenomics analysis tool; it is designed to provided maximum utility to the user by incorporating a number of analysis modules for the quantification of not only bacterial strains but also gene families and biosynthetic gene clusters. deepStrain also incorporates quality control modules and visualization tools.

---

# Citation: #
 Xinyang Zhang, Tyson Dawson, Keith A. Crandall, Ali Rahnavard (2022+), **deepStrain: jointly profile microbial strains, genes, and biosynthetic gene clusters from metagenomics data**, https://github.com/omicsEye/deepStrain

# Attention # 
----

Please check our [omicsEye Support Forum](https://forum.omicsEye.org) for common questions before open issue thread there.

----
# deepStrain user manual

## Contents ##
* [Features](#features)
* [deepStrain](#deepStrain)
    * [Requirements](#requirements)
    * [Installation](#installation)
* [Getting started with deepStrain](#getting-started-with-deepStrain)
    * [Test deepStrain](#test-deepStrain)
    * [Options](#options) 
    * [Input](#input)
    * [Output](#output)  
* [deepStrain piplines](#deepStrain-piplines)
   * [Taxonomic profiling](#taxonomic-profiling)
   * [Biosynthetic gene clusters profiling](Biosynthetic-gene-clusters-profiling)
   * [Gene families profiling](#Gene-families-profiling)
   * [Strain profiling](#Strain-profiling)
* [Utilities](#Utilities)
   * [Joint profiles](#joint-profiles)
   * [QC Visualization](#qc-visulization-demo) 
* [Real world examples](#real-world-examples)
    * [Visualization](#visulization-demo)
* [Support](#Support)
------------------------------------------------------------------------------------------------------------------------------
# Features #
1. Generality: deepStrain uses sequence reads as input with filtering and QC.

2. Mapping database
    * Taxonomic Reference Genomes
    * Gene Families
    * Biosynthesis gene clusters

3. Downstream Analysis:
    * Downstream Analysis
    * Gene Family Pathway Analysis
    * Sequence alignment & Key mutation ID

4. High Resolution Phylogenomics:
    * Identification of key species
    * Inclusion of wider diversity of strains for key species
    * Cladogram generation w/ annotations


![overall-fig](https://github.com/omicsEye/deepStrain/blob/main/img/fig1_general_pipeline_v2.png)
    
# deepStrain #

## REQUIREMENTS ##
* [Matplotlib](http://matplotlib.org/)
* [Python 3.*](https://www.python.org/download/releases/)
* [Numpy 1.9.*](http://www.numpy.org/)
* [Pandas (version >= 0.18.1)](http://pandas.pydata.org/getpandas.html)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [bowtie2 v2.4.4](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [biopython 1.79](https://biopython.org/)
* [diamond 2.0.14](https://github.com/bbuchfink/diamond)

## INSTALLATION ##
 - not finished yet

Linux based and Mac OS:
* First open a terminal 
```
$ sudo pip3 install deepStrain
```
If you use `sudo` then you need provide admin password and teh software will be installed for all users.

You can also install it as on user home directory by providing `--user` or specifying a path by providing a pATH AFTER `-t` option.

Windows OS:
* First open a Command Prompt terminal as administrator 
then run the following command 

```
$ pip3 install deepStrain
```

* You can replace `pip3` by `pip` if you have only Python 3 installed on your computer. `pip3` specifies to install `deepStrain` for Python 3. 

------------------------------------------------------------------------------------------------------------------------------

# Getting Started with deepStrain #

## Test deepStrain ##

To test if deepStrain is installed correctly, you may run the following command in the terminal:

```
#!cmd

deepStrain -h

```
Which yields deepStrain command line options.


## Options ##

```
$ deepStrain -h
usage: deepStrain [-h] [--version] [-i INPUT] -o OUTPUT 
 
optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -i INPUT, --input INPUT
                        the input file are reads with FASTA/FASTQ format
                         
  -o OUTPUT, --output OUTPUT
                        the output directory
  -v, --verbose         additional output is printed
```


## Input ##

The two required input parameters are:

1. ``-i or --input:`` reads.
2. ``--output-folder``: a folder containing all the output files

A list of all options are provided in #options section. 

## Output ##
```
$ deepStrain -h
usage: deepStrain [-h] 

```
# deepStrain piplines # 
## Taxonomic profiling ##
1. Bayesian Reassignment


![tax1](https://github.com/omicsEye/deepStrain/blob/main/img/taxProfile1.png)

2. Taxonomic profiling


![tax2](https://github.com/omicsEye/deepStrain/blob/main/img/taxProfile2.png)


3. Visualization


![tax3](https://github.com/omicsEye/deepStrain/blob/main/img/taxProfile3.png)

# Real world example #
## Visulization Demo ##
1. Go to the deepStrain/Notebooks, download FiveTargetNum.tsv, FiveTargetReads.tsv and stackedplot.ipynb.
2. FiveTargetNum.tsv and FiveTargetReads.tsv are two output files that generated from deepStrain.
3. Run the code either on the google colab or in your loacl environment.
4. The stacked bar plot show the composition distribution and their corresponding reads.
5. The final look could be liked the following:
![stacked plot](https://github.com/omicsEye/deepStrain/blob/main/Notebooks/stackedplot.png)
