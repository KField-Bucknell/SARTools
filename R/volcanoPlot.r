#' Volcano plots
#'
#' Volcano plot for each comparison: -log10(adjusted P value) vs log2(FC) with one dot per feature (red dot for a differentially expressed feature, black dot otherwise)
#'
#' @param complete A \code{list} of \code{data.frame} containing features results (from \code{exportResults.DESeq2()} or \code{exportResults.edgeR()})
#' @param alpha cut-off to apply on each adjusted p-value
#' @param outfile TRUE to export the figure in a png file
#' @param adjPlimit optional limit for -log10 adjusted P value
#' @return A file named volcanoPlot.png in the figures directory containing one volcano plot per comparison
#' @author Hugo Varet

volcanoPlot <- function(complete, alpha=0.05, outfile=TRUE, adjPlimit=0){
  ncol <- ifelse(length(complete)<=4, ceiling(sqrt(length(complete))), 3)
  nrow <- ceiling(length(complete)/ncol)
  if (outfile) png(filename="figures/volcanoPlot.png", width=cairoSizeWrapper(1800*ncol), height=cairoSizeWrapper(1800*nrow), res=300)
    par(mfrow=c(nrow,ncol))
    for (name in names(complete)){
      complete.name <- complete[[name]]
      complete.name$padj[which(complete.name$padj==0)] <- .Machine$double.xmin
      log10pval <- -log10(complete.name$padj)
      if (adjPlimit==0) {
        ylim <- c(0,1) * max(quantile(-log10(complete.name$padj), probs=0.99, na.rm=TRUE), 3)
      } else {
        ylim <- c(0,adjPlimit)
      }
      plot(complete.name$log2FoldChange, pmin(ylim[2], log10pval), ylim=ylim, las=1, cex=0.45,
           xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~adjusted~P~value),
           col=ifelse(complete.name$padj <= alpha, "red", "black"), pch=ifelse(log10pval >= ylim[2], 2, 20),
           main=paste0("Volcano plot - ",gsub("_"," ",name)))
      abline(h=-log10(alpha), lty=2, col="lightgray")
    }
  if (outfile) dev.off()
}
