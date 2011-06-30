# gseconvert

The purpose of this utility is to take raw, heterogenous gene
expression data (GSE files) from GEO and output a normalized matrix,
with genes as columns and experiments as rows, for conducting gene
expression meta-analyses such as described in some of
our papers:

*   [GAMMA: Global Microarray Meta-analysis](http://bioinformatics.oxfordjournals.org/content/early/2009/05/15/bioinformatics.btp290)

## Requirements

*  R language (tested on 2.13.0), including RScript, on your PATH
*  GNU Make and GNU grep (needs to support -P option)
*  ~150 GB free disk space

## Usage

First, download GSEs from NCBI's FTP servers,
GEO and Entrez Gene metadata, and R package dependencies by running:

    make download

Alternatively if you've already downloaded the GSEs, you can symlink
`data/GSE` to a flat directory containing the GSEs and download the
other dependencies individually (see the Makefile).

Now make the matrix for the species of your choice by running, e.g.: 

    make matrix SPECIES="Homo sapiens"

being sure to surround the species of interest by quotes.

## License

Copyright (C) 2011 Oklahoma Medical Research Foundation

Distributed under the Eclipse Public License.
