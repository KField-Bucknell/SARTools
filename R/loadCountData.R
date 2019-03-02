#' Load count files
#'
#' Load one count file per sample thanks to the file names in the target file.
#'
#' @param target target \code{data.frame} of the project returned by \code{loadTargetFile()}
#' @param rawDir path to the directory containing the count files
#' @param skip number of lines of the data file to skip before beginning to read data
#' @param idColumn the column number with feature Ids (default = 1)
#' @param countColumn the column number with counts (default = 2)
#' @param header default set to header=TRUE, set to FALSE if count file has no header row
#' @param featuresToRemove vector of feature Ids (or character string common to feature Ids) to remove from the counts
#' @return The \code{matrix} of raw counts with row names corresponding to the feature Ids and column names to the sample names as provided in the first column of the target.
#' @details If \code{featuresToRemove} is equal to \code{"rRNA"}, all the features containing the character string "rRNA" will be removed from the counts.
#' @author Marie-Agnes Dillies and Hugo Varet

loadCountData <- function(target, rawDir="raw", skip=0, header=TRUE, idColumn=1, countColumn=2, featuresToRemove=c("alignment_not_unique", "ambiguous", "no_feature", "not_aligned", "too_low_aQual")){
  
  labels <- as.character(target[,1])
  files <- as.character(target[,2])
  
  # detect if input columns have exected format
  f1 <- read.table(file.path(rawDir, files[1]), sep="\t", quote="\"", header=header, skip=skip, stringsAsFactors=FALSE)
  if (!(is.character(f1[,idColumn]) & is.numeric(f1[,countColumn]) & (ncol(f1) >= max(idColumn,countColumn)) )){
    stop("Columns do not match expectations. Confirm that column ", idColumn, 
         " contains feature Ids and column ", countColumn, " contains counts.")
  }
  
  rawCounts <- read.table(file.path(rawDir, files[1]), sep="\t", quote="\"", header=header, skip=skip, stringsAsFactors=FALSE)
  rawCounts <- rawCounts[,c(idColumn, countColumn)]
  colnames(rawCounts) <- c("Id", labels[1])
  if (any(duplicated(rawCounts$Id))) stop("Duplicated feature names in ", files[1])
  cat("Loading files:\n")
  cat(files[1],": ",length(rawCounts[,labels[1]])," rows and ",sum(rawCounts[,labels[1]]==0)," null count(s)\n",sep="")
  
  for (i in 2:length(files)){
    tmp <- read.table(file.path(rawDir, files[i]), sep="\t", quote="\"", header=header, skip=skip, stringsAsFactors=FALSE)
    tmp <- tmp[,c(idCol, countsCol)]
    colnames(tmp) <- c("Id", labels[i])
    if (any(duplicated(tmp$Id))) stop("Duplicated feature names in ", files[i])
    rawCounts <- merge(rawCounts, tmp, by="Id", all=TRUE)
    cat(files[i],": ",length(tmp[,labels[i]])," rows and ",sum(tmp[,labels[i]]==0)," null count(s)\n",sep="")
  }
  
  rawCounts[is.na(rawCounts)] <- 0
  counts <- as.matrix(rawCounts[,-1])
  rownames(counts) <- rawCounts[,1]
  counts <- counts[order(rownames(counts)),]
  
  # check that input counts are integers to fit edgeR and DESeq2 requirements
  # round them off to integers if they are not
  if (any(counts %% 1 != 0)) {
    counts <- as.integer(round(counts))
    cat("Rounded counts to integers to meet requirements of edgeR and DESeq2.\n")
  }
  
  cat("\nTop of the counts matrix:\n")
  print(head(counts))
  cat("\nBottom of the counts matrix:\n")
  print(tail(counts))
  return(counts)
}
