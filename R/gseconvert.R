library(GEOquery)
library(GEOmetadb)
library(preprocessCore)
library(stringr)
library(memoise)
library(plyr)

options(warn=-1)
ensure.geometadb <- function() {
   if (!file.exists("GEOmetadb.sqlite"))
     GEOmetadb::getSQLiteFile()
}

## Functions for mapping among IDs
query_geometadb <- function(qry) {
  conn <- dbConnect(SQLite(), "GEOmetadb.sqlite")
  sqliteQuickSQL(conn, qry) #TODO close
}

get_platforms_for_species <- function(species) {
  ##Currently restricts to one-color
  query_geometadb(sprintf("SELECT gpl,count(gpl) as count FROM gsm
     WHERE organism_ch1='%s'
     AND channel_count=1
     GROUP BY gpl
     ORDER BY count(gpl) DESC", species))$gpl
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

postprocess.matrix <- function(m) {
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

  return(round(m[qcrows,]))
}

GENE_COLS <- c("ENTREZ_GENE_ID", "GeneID", "Entrez_Gene_ID", "GENE", "Gene")

eset_to_matrix <- function(eSet) {
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
  
  #Sort columns
  m <- m[,sort(colnames(m))]
  
  #postprocess
  result <- postprocess.matrix(m)

  #Remove all rows that all purely NA
  result <- result[,apply(result,1,function(row) !all(is.na(row)))]
  return(result)
}

read.gse<- function(gseAcc, platform=NA, base_dir="data/GSE/") {
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


## Functions needed later for quantile normalization
update_quantile_normalization_vector <- function(m) {

}

GSMS_USED <- c()


write.platform.matrix <- function(platform, outfile) {
  ensure.geometadb()
  
  QNORM_COUNTS <- NULL
  QNORM_MEANS <- NULL

  gselist <- geoConvert(platform, out_type="GSE",sqlite_db_name="GEOmetadb.sqlite")[[1]]$to_acc 
  for (gseAcc in gselist) {
    gse <- read.gse(gseAcc,platform=platform)
    if (!is.na(gse)) {
      if (!is.list(gse)) {
        gse <- list(gse)
      }
      for (eSet in gse) {
          if (annotation(eSet) == platform) {
            try({
              m <- eset_to_matrix(eSet)
              m <- m[setdiff(rownames(m),GSMS_USED),] ##Don't reuse GSMs
              if (!empty(m)) {
                GSMS_USED <- c(GSMS_USED,rownames(m))
                #update the qunatile normalization vector
                if (is.null(QNORM_COUNTS)) {
                  QNORM_COUNTS <- rep(0,ncol(m))
                  QNORM_MEANS <- rep(-1,ncol(m))
                }
                for (i in 1:nrow(m)) {
                  v <- sort(m[i,]) ##Removes NAs also
                  range <- (ncol(m)-length(v)):length(v) ##Right side of the vector
                  QNORM_MEANS[range] <- ((QNORM_MEANS[range] * QNORM_COUNTS[range]) + v) / (QNORM_COUNTS[range] + 1)
                  QNORM_COUNTS[range] <- QNORM_COUNTS[range] + 1
                }
                write.table(m, file=outfile,
                            col.names=!file.exists(outfile), append=TRUE, sep="\t")
              }
            })
        }
      }
    }
  }
  return(QNORM_MEANS)
}


normalize.matrix <- function(f.in, f.out, qnorm.means) {
  n_rows <- as.numeric(strsplit(system(paste("wc -l", f.in), intern=T), " ")[[1]][1])
  
  ## Quantile normalize the big matrix, a chunk at a time
  genes <- sapply(strsplit(readLines(f.in,n=1), "\t")[[1]], function(x) substr(x,2,nchar(x)-1))
  STEP_SIZE <- 2000
  for (i in seq(1, n_rows, by=STEP_SIZE)) {
    m <- read.table(f.in,sep="\t",skip=1+STEP_SIZE*(i-1), header=F, nrow=STEP_SIZE, row.names=1)
    colnames(m) <- genes
    m <- quantile_normalize(m, qnorm.means)
    write.table(m, file=f.out, sep="\t", append=!file.exists(f.out))
  }

  ## Remove empty columns
  #system(sprintf("python src/remove_missing.py %s > %s", path_normalized, path_final))
  #unlink(path)
  #unlink(path_normalized)
  #file.rename(path_final, path)
}


##API
convert.gse.files <- function(platform, outfile) {
  tmp <- tempfile()
  qnorm.means <- write.platform.matrix(platform, tmp)
  normalize.matrix(tmp, outfile, qnorm.means)
}
#write.platform.matrix(platform, outfile) <- 
#write.species.matrix(platform, outfile)
#-postprocess.matrix
