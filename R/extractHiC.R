#' @export
#sourceCpp('~/Desktop/straw-master/R/straw-R.cpp')

extractHiC <- function(hic, chrom, chromstart, chromend, resolution, plottype = "square", zrange=NULL, norm="NONE", altchrom = NULL, altchromstart = NULL, altchromend = NULL){
  ## Parse region input for straw
  
  
  #parameter = region
  # regionSplit <- unlist(strsplit(unlist(strsplit(region, ":")), "-"))
  # print(regionSplit)
  # regionCleaned <- gsub(pattern = "chr|chrom|CHR|CHROM", replacement = "", x = regionSplit)
  # print(regionCleaned)
  # regionStraw <- paste(regionCleaned, collapse = ":")
  
  regionChrom <- gsub(pattern = "chr", replacement = "", x = chrom)
  regionStraw <- paste(regionChrom, chromstart, chromend, sep = ":")

  ## Extract upper triangular using straw; reverse for lower triangular; combine unique
  if(is.null(altchrom)){
    upper <- straw_R(sprintf("%s %s %s %s BP %i", norm, hic, regionStraw, regionStraw, resolution))

  } else {
    
    regionChrom2 <- gsub(pattern = "chr", replacemtn = "", x = altchrom)
    regionStraw2 <- paste(regionChrom2, altchromstart, altchromend, sep = ":" )
    upper <- straw_R(sprintf("%s %s %s %s BP %i", norm, hic, regionStraw, regionStraw2, resolution))
    
  }
  
  ## Complete missing values; convert missing values to 0 counts
  if(plottype == "square"){
    lower <- upper[, c(2,1,3)]
    colnames(lower) <- c("x", "y", "counts")
    combined <- unique(rbind(upper, lower))
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
    
    return(as.data.frame(combinedComplete))
    
  } else {
  
    ## Scale lower and upper bounds using zrange
    if(is.null(zrange)){
      zrange <- c(0, max(upper$counts))
      upper$counts[upper$counts <= zrange[1]] <- zrange[1]
      upper$counts[upper$counts >= zrange[2]] <- zrange[2]
    } else {
      stopifnot(is.vector(zrange), length(zrange) == 2, zrange[2] > zrange[1])
      upper$counts[upper$counts <= zrange[1]] <- zrange[1]
      upper$counts[upper$counts >= zrange[2]] <- zrange[2]

    }
    
    return(as.data.frame(upper))
  }


}


#hicdata <- extractHiC("/Users/phanstiel7/Documents/mega_inter_30.hic", chrom = "chr8", chromstart = chromstart, chromend = chromend, resolution = 5000, zrange=c(0,12), plottype="triangle")


