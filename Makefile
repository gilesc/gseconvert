species-matrix : 
	Rscript src/gseconvert.R "$(SPECIES)"

platform-matrix : 
	Rscript src/gseconvert.R "$(PLATFORM)"

download-gse :
	mkdir -p data/GSE
	cd data/GSE && wget -nd -r ftp://ftp.ncbi.nih.gov/pub/geo/DATA/

download-geometadb : 
	cd data && wget http://gbnci.abcc.ncifcrf.gov/geo/GEOmetadb.sqlite.gz && gunzip GEOmetadb.sqlite.gz

download-taxonomy : 
	mkdir -p data/taxonomy
	cd data/taxonomy && wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && \
		gunzip taxdump.tar.gz && \
		tar xf taxdump.tar && \
		ls | grep -v names | xargs rm && \
		mv names.dmp ../taxonomy.dat && \
		cd .. && rm -r taxonomy

download-gene-info :
	cd data && wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz && \
		gunzip gene_info.gz

download-R-deps :
	Rscript src/installdeps.R

download : download-gse download-taxonomy download-gene-info download-geometadb download-R-deps


