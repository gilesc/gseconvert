# gseconvert

The purpose of this utility is to take raw, heterogenous gene
expression data (GSE files) from GEO and output a normalized matrix,
with genes as columns and experiments as rows, for conducting gene
expression meta-analyses such as described in some of
our papers:

*   [GAMMA: Global Microarray Meta-analysis](http://bioinformatics.oxfordjournals.org/content/early/2009/05/15/bioinformatics.btp290)

## Usage

First, download GSEs and GPL annotation data from NCBI's FTP servers,
as well as GEO metadata by running:

   `make download`

Alternatively if you've already downloaded some of these, you can symlink
data/GSE or data/GPL to flat directories containing GSEs or GPLs,
respectively (they can be gzipped).

Now make the matrix of your choice by running: 

    `make $SPECIES`

where $SPECIES is "Homo sapiens", "Mus musculus", etc., surrounded by
quotes. 


## License

Copyright (C) 2011 Oklahoma Medical Research Foundation

Distributed under the Eclipse Public License.
