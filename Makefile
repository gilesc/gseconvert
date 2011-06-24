matrix : 
	@echo "TODO"	

download-gse :
	mkdir -p data/GSE
	cd data/GSE && wget -nd -r ftp://ftp.ncbi.nih.gov/pub/geo/DATA/

download-gpl :
	mkdir -p data/GPL
	cd data/GPL && wget -nd -r ftp://ftp.ncbi.nih.gov/pub/geo/DATA/annotation/platforms/

download : download-gse download-gpl


