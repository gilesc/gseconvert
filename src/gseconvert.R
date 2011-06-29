library(GEOquery)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Ce.eg.db)
library(org.EcK12.eg.db)
library(org.Rn.eg.db)
library(GEOmetadb)
library(preprocessCore)
library(impute)

postprocess.matrix <- function(m, imputeMissing=TRUE) {
  if (imputeMissing) {
    m <- impute.knn(m)$data
  }
  #As per Mikhail's paper:
  #floor the top 0.1% of genes (todo. do per experiment)
  high.genes <- names(sort(rowMeans(m, na.rm=TRUE), decreasing=TRUE))
  high.genes <- high.genes[1:ceiling(length(high.genes) * 0.001)]
  m[high.genes,] <- apply(m[high.genes,], 2, min)
  #scale to 10000 by experiment
  m <- sweep(m, MARGIN=2, 10000.0 / apply(m, 2, function(col) max(col, na.rm=TRUE)), `*`)
  #quantile normaliztion
  axes <- dimnames(m)
  m <- normalize.quantiles(m, copy=TRUE)
  dimnames(m) <- axes
  return(round(m))
}


geneCols <- c("ENTREZ_GENE_ID", "GeneID", "Entrez_Gene_ID", "GENE", "Gene")

eSet2Matrix <- function(eSet, allGenes) {
  pd <- pData(featureData(eSet))
  col <- geneCols[geneCols %in% names(pd)][1]
  genes <- as.vector(pd[,col])
  m <- t(exprs(eSet))
  colnames(m) <- genes
  m <- m[,colnames(m) != ""]
  
  #Handle columns with multiple Entrez IDs
  cols <- grep("///", colnames(m))
  subMs <- sapply(colnames(m)[cols], 
            function(cn) { 
              cnt <- length(gregexpr("///",cn)[[1]]) + 1
              matrix(rep(m[,cn], cnt), ncol=cnt)
            })
  subM <- NULL
  for (sM in subMs) {
    subM <- cbind(subM, sM)
  }
  colnames(subM) <- unlist(strsplit(colnames(m)[cols], " /// "))
  
  
  #Combine normal and "abnormal" columns
  m <- cbind(m[,!(seq(1,ncol(m)) %in% cols)], subM)
  
  #Take most active probe (simplistic approach needs to be changed in the future..)
  m <- m[,order(-colMeans(m))]
  m <- m[,sort(unique(colnames(m)))]
  
  #Put this new matrix into a bigger matrix containing a column for all genes in the species
  m <- m[,intersect(colnames(m),allGenes)]
  result <- matrix(nrow=nrow(m),ncol=length(allGenes))
  rownames(result) <- rownames(m)
  colnames(result) <- allGenes
  idxs <- match(colnames(m), allGenes)
  result[,idxs] <- m
  colnames(result) <- paste("LL:", colnames(result), sep="")
  return(result)
}


base_dir <- "data/GSE/"
readGSE <- function(gseAcc, platform=NA) {
  gse <- NA
  path <- ""
  if (!is.na(platform)) {
    path <- paste(base_dir, gseAcc, "-", platform, "_series_matrix.txt.gz", sep="")
  }
  if (!file.exists(path)) {
    path <- paste(base_dir, gseAcc, "_series_matrix.txt.gz", sep="")
  }
  if (file.exists(path)) {
    try(
      gse <- getGEO(filename=path))
  } else {
    try(
      gse <- getGEO(gseAcc))
  }
  return(gse)
}

result <- NA

writePlatformMatrix <- function(platform) {
  outfile <- paste("/home/gilesc/Desktop/",platform, ".tsv", sep="")
  genes <- getAllGenesForSpecies##TODOTODOTODOTODODOTODODODODODO
  gselist <- geoConvert(platform, out_type="GSE")[[1]]$to_acc 
  for (gseAcc in gselist) {
    print(gseAcc)
    gse <- readGSE(gseAcc,platform=platform)
    if (!is.na(gse)) {
      if (!is.list(gse)) {
      gse <- list(gse)
      }
      for (eSet in gse) {
          if (annotation(eSet) == platform) {
            try({
              m <- eSet2Matrix(eSet, genes)
              write.table(m, file=outfile, col.names=!file.exists(outfile), append=TRUE, sep="\t")})
        }
      }
    }
  }
  return(result)
}

writeMatrix <- function(species) {
  ##TODO:
}


##accessions <- c("GPL96", "GPL570", "GPL97", "GPL571") #human
##accessions <- c("GPL1261", "GPL81", "GPL339", "GPL3667", "GPL8321", "GPL6885") #mouse
##accessions <- c("GPL200") #c elegans
##accessions <- c("GPL199", "GPL3154") #e. coli
accessions <- c("GPL1355", "GPL85", "GPL341", "GPL4135", "GPL6101") #rat
for (acc in accessions) {
  writePlatformMatrix(acc) 
}

