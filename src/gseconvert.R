library(GEOquery)
library(GEOmetadb)
library(preprocessCore)
library(stringr)
library(memoise)
library(plyr)

options(warn=-1)


## Functions for mapping among IDs
query_geometadb <- function(qry) {
    conn <- dbConnect(SQLite(), "data/GEOmetadb.sqlite")
    sqliteQuickSQL(conn, qry) #TODO close
}

get_species_for_platform <- function(platform) {
  query_geometadb(sprintf("SELECT organism from gpl
    WHERE gpl='%s' LIMIT 1", platform))$organism
}

get_platforms_for_species <- function(species) {
  ##Currently restricts to one-color
  query_geometadb(sprintf("SELECT gpl,count(gpl) as count FROM gsm
     WHERE organism_ch1='%s'
     AND channel_count=1
     GROUP BY gpl
     ORDER BY count(gpl) DESC", species))$gpl
}


## Functions needed later for quantile normalization
N_ROWS_REJECTED_QC <- 0

QNORM_COUNTS <- NULL
QNORM_MEANS <- NULL
update_quantile_normalization_vector <- function(m) {
  if (is.null(QNORM_COUNTS)) {
    QNORM_COUNTS <<- rep(0,ncol(m))
    QNORM_MEANS <<- rep(-1,ncol(m))
  }
  for (i in 1:nrow(m)) {
    v <- sort(m[i,]) ##Removes NAs also
    range <- (ncol(m)-length(v)):length(v) ##Right side of the vector
    QNORM_MEANS[range] <<- ((QNORM_MEANS[range] * QNORM_COUNTS[range]) + v) / (QNORM_COUNTS[range] + 1)
    QNORM_COUNTS[range] <<- QNORM_COUNTS[range] + 1
  }
}

quantile_normalize <- function(m,v) {
  ##Where v is the quantile normalization mean vector (QNORM_MEANS)
  v[v==-1] <- NA
  t(apply(m,1,function(row) {
    row <- sort(row,na.last=FALSE)
    result <- v
    names(result) <- names(row)
    result[colnames(m)]
  }))
}

postprocess_matrix <- function(m) {
  m <- t(apply(m,1, function(row) {
    min.val <- min(row,na.rm=T)
    max.val <- max(row,na.rm=T)
    ## Identify experiments that are log-based and rescue them 
    if (min.val >= 7 & max.val <= 16) {
      row <- 2 ^ row
      min.val <- min(row,na.rm=T)
      max.val <- max(row,na.rm=T)
    }
                 
    ##scale to 10000 by experiment
    10000 * (row - min.val) / max(1, (max.val - min.val))
  }))
  
  #As per Mikhail's paper:
  #Floor the top 0.1% of genes
  high.genes <- names(sort(colMeans(m, na.rm=TRUE), decreasing=TRUE))
  high.genes <- high.genes[1:ceiling(length(high.genes) * 0.001)]
  m[,high.genes] <- apply(m[,high.genes], 1, min)
  
  ## Quality control criteria:
  ## 1. Mean and median >= 0
  ## 2. Mean-median ratio >= 1.2
  ## 3. Any negative values
  qcrows <- apply(m,1,function(row) {
    row <- row[!is.na(row)]
    row.mean <- mean(row)
    row.median <- median(row)
    (row.median >= 0) & (row.mean >= 0) & (row.mean / row.median >= 1.2) & (all(row>=0)) & (length(row) > 0)
  })

  N_ROWS_REJECTED_QC <<- N_ROWS_REJECTED_QC + sum(!qcrows,na.rm=T)
  return(round(m[qcrows,]))
}

GENE_COLS <- c("ENTREZ_GENE_ID", "GeneID", "Entrez_Gene_ID", "GENE", "Gene")

eset_to_matrix <- function(eSet, allGenes) {
  pd <- pData(featureData(eSet))
  col <- GENE_COLS[GENE_COLS %in% names(pd)][1]
  if (empty(pd) || is.null(col)) {
    ##This platform is not supported
    print("WARNING: This platform is not supported! Recovering...")
    return(c())
  }
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
  
  #Take most active probe
  m <- m[,order(-colMeans(m))]
  m <- m[,sort(unique(colnames(m)))]
  
  #Put this new matrix into a bigger matrix containing a column for all genes in the species
  m <- m[,intersect(colnames(m),allGenes)]
  result <- matrix(nrow=nrow(m),ncol=length(allGenes))
  rownames(result) <- rownames(m)
  colnames(result) <- allGenes
  idxs <- match(colnames(m), allGenes)
  result[,idxs] <- m

  #postprocess
  result <- postprocess_matrix(result)

  #Remove all rows that all purely NA
  result <- result[,apply(result,1,function(row) !all(is.na(row)))]
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

get_name <- function(species_or_platform) {
  #Get "subname" of file
  name <- ifelse(substr(species_or_platform,1,3)=="GPL", species_or_platform,
         tolower(species_or_platform))
  str_replace_all(name, " ", "_")
}
get_outfile <- function(species_or_platform, suffix=NA) {
  sprintf("data/gse-%s%s.mtx",
          get_name(species_or_platform),
          ifelse(is.na(suffix), "", sprintf("-%s", suffix)))
}

GSMS_USED <- c()
append_platform_to_matrix_file <- function(platform, outfile=get_outfile(platform)) {
  genes <- get_genes_for_species(get_species_for_platform(platform))
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
              m <- m[setdiff(rownames(m),GSMS_USED),] ##Don't reuse GSMs
              if (!empty(m)) {
                GSMS_USED <- c(GSMS_USED,rownames(m))
                update_quantile_normalization_vector(m)
                write.table(m, file=outfile,
                            col.names=!file.exists(outfile), append=TRUE, sep="\t")
              }
            })
        }
      }
    }
  }
}

normalize <- function(platform_or_species) {
  cat(paste(N_ROWS_REJECTED_QC, "experiments (GSMs) were rejected due to failing QC.\n"))
  path <- get_outfile(platform_or_species)
  n_rows <- as.numeric(strsplit(system(paste("wc -l", path), intern=T), " ")[[1]][1])
  
  ## Quantile normalize the big matrix, a chunk at a time
  genes <- sapply(strsplit(readLines(path,n=1), "\t")[[1]], function(x) substr(x,2,nchar(x)-1))
  STEP_SIZE <- 2000
  for (i in seq(1, n_rows, by=STEP_SIZE)) {
    m <- read.table(path,sep="\t",skip=1+STEP_SIZE*(i-1), header=F, nrow=STEP_SIZE, row.names=1)
    colnames(m) <- genes
    m <- quantile_normalize(m, QNORM_MEANS)
    path_normalized <- get_outfile(platform_or_species,suffix="normalized")
    path_final <- get_outfile(platform_or_species, suffix="final")
    write.table(m, file=path_normalized, sep="\t", append=!file.exists(path_normalized))
  }

  ## Remove empty columns
  system(sprintf("python src/remove_missing.py %s > %s", path_normalized, path_final))
  unlink(path)
  unlink(path_normalized)
  file.rename(path_final, path)
}

write_platform_matrix <- function(platform) {
  outfile <- get_outfile(platform)
  file.remove(outfile)
  append_platform_to_matrix_file(platform)
  normalize(platform)
  
}

write_species_matrix <- function(species) {
  outfile <- get_outfile(species)
  file.remove(outfile)
  for (platform in get_platforms_for_species(species)[1:10]) {
    ##For now, the rarer platforms seem to be more trouble than they're worth, so using only the top 10
    print(platform)
    append_platform_to_matrix_file(platform,outfile=get_outfile(species))
  }
  normalize(species)
}

args <- as.character(commandArgs(trailingOnly=T))
if (tolower(substr(args,1,3))=="gpl") {
  write_platform_matrix(args)
} else {
  write_species_matrix(args)
}

