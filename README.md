# gseconvert

The purpose of this utility is to take raw, heterogenous gene
expression data (GSE files) from GEO and output a normalized matrix,
with genes as columns and experiments as rows, for conducting gene
expression meta-analyses such as described in some of
our papers:

*   [GAMMA: Global Microarray Meta-analysis](http://bioinformatics.oxfordjournals.org/content/early/2009/05/15/bioinformatics.btp290)

## Usage

First, download GSEs from NCBI's FTP servers,
as well as GEO metadata by running:

    make download

Alternatively if you've already downloaded some of these, you can symlink
`data/GSE` a flat directory containing GSEs, (they can be gzipped).

Now make the matrix for the species of your choice by running, e.g.: 

    make matrix SPECIES="Homo sapiens"

being sure to surround the species of interest by quotes.

## License

Copyright (C) 2011 Oklahoma Medical Research Foundation

Distributed under the Eclipse Public License.
