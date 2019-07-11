#' @export

extractHiC <- function(hic, region, resolution, zrange=NULL, norm="NONE"){
  ## Parse region input for straw
  regionSplit <- unlist(strsplit(unlist(strsplit(region, ":")), "-"))
  regionCleaned <- gsub(pattern = "chr|chrom|CHR|CHROM", replacement = "", x = regionSplit)
  regionStraw <- paste(regionCleaned, collapse = ":")
  
  ## Extract upper triangular using straw; reverse for lower triangular; combine unique
  upper <- straw_R(sprintf("%s %s %s %s BP %i", norm, hic, regionStraw, regionStraw, resolution))
  lower <- upper[, c(2,1,3)]
  colnames(lower) <- c("x", "y", "counts")
  combined <- unique(rbind(upper, lower))
  
  ## Complete missing values; convert missing values to 0 counts
  combinedComplete <- tidyr::complete(combined, x, y)
  combinedComplete$counts[is.na(combinedComplete$counts)] <- 0
  
  ## Scale lower and upper bounds using zrange
  if(is.null(zrange)){
    zrange <- c(0, max(combinedComplete$counts))
    combinedComplete$counts[combinedComplete$counts <= zrange[1]] <- zrange[1]
    combinedComplete$counts[combinedComplete$counts >= zrange[2]] <- zrange[2]
  } else {
    stopifnot(is.vector(zrange), length(zrange) == 2, zrange[2] > zrange[1])
    combinedComplete$counts[combinedComplete$counts <= zrange[1]] <- zrange[1]
    combinedComplete$counts[combinedComplete$counts >= zrange[2]] <- zrange[2]
  }
  
  return(combinedComplete)
}
