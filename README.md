# seqSight #

**seqSight** is a tool to jointly profile microbial strains, genes, and biosynthetic gene clusters from metagenomics data; it is designed to provided maximum utility to the user by incorporating a number of analysis modules for the quantification of not only bacterial strains but also gene families and biosynthetic gene clusters. seqSight also incorporates quality control modules and visualization tools.
 
---

# Citation: #
 Xinyang Zhang, Tyson Dawson, Keith A. Crandall, Ali Rahnavard (2022+), **seqSight: jointly profile microbial strains, genes, and biosynthetic gene clusters from metagenomics data**, https://github.com/omicsEye/seqSight

# Attention # 
----

Please check our [omicsEye Support Forum](https://forum.omicsEye.org) for common questions before open issue thread there.

----
# seqSight user manual

## Contents ##
* [Features](#features)
* [seqSight](#seqSight)
    * [Requirements](#requirements)
    * [Installation](#installation)
* [Getting started with seqSight](#getting-started-with-seqSight)
    * [Test seqSight](#test-seqSight)
    * [Options](#options) 
    * [Input](#input)
    * [Output](#output)  
* [seqSight pipelines](#seqSight-pipelines)
   * [Taxonomic profiling](#taxonomic-profiling)
   * [Biosynthetic gene clusters profiling](Biosynthetic-gene-clusters-profiling)
* [Utilities](#Utilities)
   * [Join tables](#join-tables)
   * [QC Visualization](#qc-visulization-demo)
* [BGC Database](#BGC-Database)
 
| Database | Description |
| --- | --- |
| `MIBiG` | MIBiG is a comprehensive database that focuses on natural product BGCs. It focuses on providing curated BGCs with associated metadata, including chemical structures, gene annotations, and experimental data. MIBiG serves as a centralized resource for **known BGCs** and aims to standardize the annotation and reporting of BGC information..|
| `antiSMASH` | antiSMASH, in addition to being a BGC prediction tool, maintains a database of **predicted BGCs**. The database includes a broad range of BGCs, covering various secondary metabolite classes.|
| `ClusterFinder` | The ClusterFinder database is a repository of **predicted BGCs** identified by the ClusterFinder tool. It focuses on BGCs identified in **bacterial genomes** and provides information on gene clusters, predicted products, and **associated metadata**. It aims to facilitate the exploration and analysis of BGC diversity in **bacteria**.|
| `PRISM` | PRISM hosts a database of **predicted BGCs** identified by the PRISM tool. It includes a collection of BGCs and associated secondary metabolite annotations. The database aims to facilitate the analysis and comparison of BGCs **across different genomes**.|
| `BiG-FAM` | BiG-FAM maintains a database of biosynthetic gene families identified by the BiG-FAM tool. It provides information on conserved gene families across BGCs and supports comparative analysis and functional characterization.|

* [Real world examples](#real-world-examples)
    * [Visualization](#visulization-demo)
* [Support](#Support)
------------------------------------------------------------------------------------------------------------------------------
# Features #
1. Generality: seqSight uses sequence reads as input with filtering and QC.

2. Mapping database
    * Taxonomic Reference Genomes
    * Gene Families
    * Biosynthesis gene clusters

3. Downstream Analysis:
    * Downstream Analysis
    * Gene Family Pathway Analysis
    * Sequence alignment & Key mutation ID


![overall-fig](https://github.com/omicsEye/seqSight/blob/main/img/fig1_general_pipeline.png)
    
# seqSight #

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

Linux-based and Mac OS:
* First open a terminal 
```
$ sudo pip3 install seqSight
```
If you use `sudo` then you need to provide an admin password and the software will be installed for all users.

You can also install it on the user's home directory by providing `--user` or specifying a path by providing a pATH AFTER `-t` option.

Windows OS:
* First open a Command Prompt terminal as administrator 
then run the following command 

```
$ pip3 install seqSight
```

* You can replace `pip3` by `pip` if you have only Python 3 installed on your computer. `pip3` specifies to install `seqSight` for Python 3. 

<span style="color:#033C5A">*If you have a working conda on your system, you can safely skip to step three*</span>.

* Install *conda*  
Go to the [Anaconda website](https://www.anaconda.com/) and download the latest version for your operating system.  
*DO NOT FORGET TO ADD CONDA TO your system PATH*
* Second is to check for conda availability  
open a terminal (or command line for Windows users) and run:
```
conda --version
```
it should output something like this:
```
conda 4.9.2
```
<span style="color:#fc0335">if not, you must make *conda* available to your system for further steps.</span>
if you have problems adding conda to PATH, you can find instructions [here](https://docs.anaconda.com/anaconda/user-guide/faq/).
  
* Third create a new conda environment (let's call it seqSight_env) with the following command:
```
conda create --name seqSight_env python=3.8
```
* Then activate your conda environment:
```commandline
conda activate seqSight_env 
```
* Finally, install *seqSight*:

* You can directly install if from GitHub:
```commandline
python -m pip install git+https://github.com/omicsEye/seqSight
```
* or before running the following line you should change your directory to the same directory that you have cloned the 
  seqSight repo:
```commandline
python -m pip install .
```

------------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------------

# Getting Started with seqSight #

## Test seqSight ##

To test if seqSight is installed correctly, you may run the following command in the terminal:

```
#!cmd

seqSight -h

```
Which yields seqSight command line options.


## Options ##

```
$ seqSight -h
usage: seqSight [-h] [--version] [-i INPUT] -o OUTPUT 
 
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
$ seqSight -h
usage: seqSight [-h] 

```
# seqSight piplines # 
## Taxonomic profiling ##
1. Bayesian Reassignment


![tax1](https://github.com/omicsEye/seqSight/blob/main/img/taxProfile1.png)

2. Taxonomic profiling


![tax2](https://github.com/omicsEye/seqSight/blob/main/img/taxProfile2.png)


3. Visualization


![tax3](https://github.com/omicsEye/seqSight/blob/main/img/taxProfile3.png)

## Utilities ##

### Join tables ###



## Visulization Demo ##
1. Go to the seqSight/Notebooks, download FiveTargetNum.tsv, FiveTargetReads.tsv and stackedplot.ipynb.
2. FiveTargetNum.tsv and FiveTargetReads.tsv are two output files that generated from seqSight.
3. Run the code either on the google colab or in your loacl environment.
4. The stacked bar plot show the composition distribution and their corresponding reads.
5. The final look could be liked the following:
![stacked plot](https://github.com/omicsEye/seqSight/blob/main/Notebooks/stackedplot.png)
