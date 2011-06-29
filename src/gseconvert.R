library(GEOquery)
library(GEOmetadb)
library(preprocessCore)
library(stringr)
library(memoise)

options(warn=-1)

postprocess_matrix <- function(m) {
  ## TODO: fix log experiments
  
  #As per Mikhail's paper:
  #floor the top 0.1% of genes (todo. do per experiment)
  high.genes <- names(sort(colMeans(m, na.rm=TRUE), decreasing=TRUE))
  high.genes <- high.genes[1:ceiling(length(high.genes) * 0.001)]
  m[,high.genes] <- apply(m[,high.genes], 1, min)
  #scale to 10000 by experiment
  m <- t(apply(m,1, function(row) {
    min.val <- min(row,na.rm=T)
    max.val <- max(row,na.rm=T)
    10000 * (row - min.val) / max(1, (max.val - min.val))
  }))
  ## Quality control criteria:
  ## 1. Mean and median >= 0
  ## 2. Mean-median ratio >= 1.2
  ## 3. <= 1% negative values
  qcrows <- apply(m,1,function(row) {
    row.mean <- mean(row,na.rm=T)
    row.median <- median(row,na.rm=T)
    (row.median >= 0) & (row.mean >= 0) & (row.mean / row.median >= 1.2) & (sum(row<0) <= length(row) / 100)
  })
  return(round(m[qcrows,]))
}

geneCols <- c("ENTREZ_GENE_ID", "GeneID", "Entrez_Gene_ID", "GENE", "Gene")

eset_to_matrix <- function(eSet, allGenes) {
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

  #postprocess
  result <- postprocess_matrix(result)
  return(result)
}

read_gse<- function(gseAcc, platform=NA, base_dir="data/GSE/") {
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
      gse <- getGEO(filename=path, AnnotGPL=T))
  } else {
    try(
      gse <- getGEO(gseAcc, AnnotGPL=T))
  }
  return(gse)
}

get_genes_for_species <- memoize(function(species) {
  tax_id <- system(sprintf("grep -P \"\t%s\t\" data/taxonomy.dat | cut -f1", species), intern=T)[1]
  system(sprintf("grep -P \"^%s\t\" data/gene_info | cut -f2", tax_id), intern=T)
})

get_outfile <- function(species) {
   sprintf("data/gse-%s.mtx",
                     str_replace_all(species, " ", "_"))
}

append_platform_to_matrix_file <- function(species, platform) {
  outfile <- get_outfile(species)

  genes <- get_genes_for_species(species)
  gselist <- geoConvert(platform, out_type="GSE",sqlite_db_name="data/GEOmetadb.sqlite")[[1]]$to_acc 
  for (gseAcc in gselist) {
    print(gseAcc)
    gse <- read_gse(gseAcc,platform=platform)
    if (!is.na(gse)) {
      if (!is.list(gse)) {
        gse <- list(gse)
      }
      for (eSet in gse) {
          if (annotation(eSet) == platform) {
            try({
              m <- eset_to_matrix(eSet, genes)
              outfile
                write.table(m, file=outfile,
                            col.names=!file.exists(outfile), append=TRUE, sep="\t")
            })
        }
      }
    }
  }
}

get_platforms_for_species <- function(species) {
  ##Currently restricts to one-color
  conn <- dbConnect(SQLite(), "data/GEOmetadb.sqlite")
  qry <- sprintf("SELECT gpl,count(gpl) as count FROM gsm
 WHERE organism_ch1='%s'
 AND channel_count=1
 GROUP BY gpl
 ORDER BY count(gpl) DESC", species)
  sqliteQuickSQL(conn, qry)$gpl
}

write_matrix <- function(species) {
  file.remove(get_outfile(species))
  for (platform in get_platforms_for_species(species)) {
    print(platform)
    append_platform_to_matrix_file(species, platform)
  }
}

species <- commandArgs(trailingOnly=T)
write_matrix(species)


