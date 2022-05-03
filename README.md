# deepStrain: #

**deepStrain** is a strain-level metagenomics analysis tool; it is designed to provided maximum utility to the user by incorporating a number of analysis modules for the quantification of not only bacterial strains but also gene families and biosynthetic gene clusters. deepStrain also incorporates quality control modules and visualization tools.

---

**Citation:**

----
# deepStrain user manual

## Contents ##
* [Features](#features)
* [deepStrain](#deepStrain)
    * [deepStrain approach](#deepStrain-approach)
    * [Requirements](#requirements)
    * [Installation](#installation)
* [Getting Started with deepStrain](#getting-started-with-deepStrain)
    * [Test deepStrain](#test-deepStrain)
    * [Options](#options) 
    * [Input](#input)
    * [Output](#output)  
* [Real world examples](#real-world-examples)
    * [Microbial species communities](#microbial-species-communities)
    * [Microbial strains](#microbial-strains)
    * [Cell line gene expressions](#cell-line-gene-expressions)
* [Support](#Support)
------------------------------------------------------------------------------------------------------------------------------
# Features #
1. Generality: deepStrain uses distance matrix as input, to allow users decide about appropriate distance metric for 
their data.

2. A simple user interface (single command driven flow)
    * The user only needs to provide a distance matrix file and a metadata file (optional)

3. A complete report including main outputs:
    * A text file of clusters and related information is provided as output in a tab-delimited file, `clusters.txt`
    * Ordination plots (PCoA, PCA, MDS, and t-SNE), heatmap,and network plot are provides for ease of interpretation
    * Discretized metadata that has been used for enrichment score calculation 
    
# deepStrain #
## deepStrain approach ##


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
usage: deepStrain [-h] [--version] [-i INPUT] -o OUTPUT [-m SIMILARITY] [--metadata METADATA] [-n ESTIMATED_NUMBER_OF_CLUSTERS] [--size-to-plot SIZE_TO_PLOT]
                [-c {single,average,complete,weighted,centroid,median,ward}] [--plot] [--resolution {high,medium,low}] [--enrichment {nmi,freq}] [-v]

Multi-resolution clustering using hierarchical clustering and Silhouette score.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -i INPUT, --input INPUT
                        the input file D*N, Rows: D features and columns: N samples OR 
                        a distance matrix file D*D (rows and columns should be the same and in the same order) 
                         
  -o OUTPUT, --output OUTPUT
                        the output directory
  -m SIMILARITY, --similarity SIMILARITY
                        similarity measurement {default spearman, options: spearman, nmi, ami, dmic, mic, pearson, dcor}
  --metadata METADATA   Rows are features and each column is a metadata
  -n ESTIMATED_NUMBER_OF_CLUSTERS, --estimated_number_of_clusters ESTIMATED_NUMBER_OF_CLUSTERS
                        estimated number of clusters
  --size-to-plot SIZE_TO_PLOT
                        Minimum size of cluster to be plotted
  -c {single,average,complete,weighted,centroid,median,ward}, --linkage_method {single,average,complete,weighted,centroid,median,ward}
                        linkage clustering method method {default = complete, options average, complete
  --plot                dendrogram plus heatmap
  --resolution {high,medium,low}
                        Resolution c .         Low resolution is good when clusters are well separated clusters.
  --enrichment {nmi,freq}
                        enrichment method.
  -v, --verbose         additional output is printed
```


## Input ##

The two required input parameters are:

1. ``-i or --input:`` reads.
Th input is a  symmetric distance matrix in a format of a tab-delimited text file of `n * n` where `n` is number of features 
(e.g. metabolites, stains, microbial species, individuals).
2. ``--output-folder``: a folder containing all the output files

Also, user can specify a metadata input to find enrichment score for each metadata 
* ``--metadata``: a tab-delimited text file with `n` rows for features names and `m` columns for metadata

A list of all options are provided in #options section. 

## Output ##



```
$ deepStrain -h
usage: omeClustviz [-h] [--metadata METADATA] [--shapeby SHAPEBY] -o OUTPUT [--size-to-plot SIZE_TO_PLOT] [--fig-size FIG_SIZE FIG_SIZE] [--point-size POINT_SIZE] [--show] adist clusters

deepStrain visualization script.

positional arguments:
  adist                 the input file D*N, Rows: D features and columns: N samples OR 
                        a distance matrix file D*D (rows and columns should be the same and in the same order) 
                         
  clusters              the input file D*N, Rows: D features and columns: N samples OR 
                        a distance matrix file D*D (rows and columns should be the same and in the same order) 
                         

optional arguments:
  -h, --help            show this help message and exit
  --metadata METADATA   metadata
  --shapeby SHAPEBY     the input file D*N, Rows: D features and columns: N samples OR 
                        a distance matrix file D*D (rows and columns should be the same and in the same order) 
                         
  -o OUTPUT, --output OUTPUT
                        the output directory
  --size-to-plot SIZE_TO_PLOT
                        Minimum size of cluster to be plotted
  --fig-size FIG_SIZE FIG_SIZE
                        width and height of plots
  --point-size POINT_SIZE
                        width and height of plots
  --show                show ordination plot before save

```
### Distance using genomics variation ###

# Real world example #



