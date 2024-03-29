#' extracts HiC data from .hic file using Straw
#'
#'
#' @param hic path to .hic file
#' @param format format of data wanted, which will depend on what kind of plot is ultimately desired; options are "sparse" and "full"
#' @param chrom if not alternative chromosome, chromosome of desired region
#' @param chromstart chromosome start position of chrom, in bp
#' @param chromend chromosome end position of chrom, in bp
#' @param resolution the width in bp of each pixel
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param norm hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#' @param resscale scale of normalization; options are "BP" and "FRAG"
#' @param altchrom if looking at region between two different chromosomes, this is the specified alternative chromsome
#' @param altchromstart if looking at region between two different chromosomes, start position of altchrom
#' @param altchromend if looking at region between two different chromsomes, end position of altchrom
#'
#'
#' @export

extractHiC <- function(hic, format, chrom, chromstart = NULL, chromend = NULL, resolution, zrange = NULL, norm = "NONE", resscale = "BP", altchrom = NULL, altchromstart = NULL, altchromend = NULL){

  # Parse chromosome and region in format for Straw
  if ((is.null(chromstart) & !is.null(chromend)) | (is.null(chromend) & !is.null(chromstart))){
    stop("Cannot have one \'NULL\' chromstart or chromend.")
  } else if (is.null(chromstart) & is.null(chromend)){
    regionStraw <- gsub(pattern = "chr", replacement = "", x = chrom)
  } else {
    regionChrom <- gsub(pattern = "chr", replacement = "", x = chrom)
    regionStraw <- paste(regionChrom, chromstart, chromend, sep = ":")
  }

  # Extract upper triangular using straw, depending on one chromsome interaction or multiple chromosome interactions
  if(is.null(altchrom)){
    upper <- straw_R(sprintf("%s %s %s %s %s %i", norm, hic, regionStraw, regionStraw, resscale, resolution))
  } else {
    if ((is.null(altchromstart) & !is.null(altchromend)) | (is.null(altchromend) & !is.null(altchromstart))){
      stop("Cannot have one \'NULL\' altchromstart or altchromend.")
    } else if (is.null(altchromstart) & is.null(altchromend)){
      regionStraw2 <- gsub(pattern = "chr", replacement = "", x = altchrom)
    } else {
      regionChrom2 <- gsub(pattern = "chr", replacement = "", x = altchrom)
      regionStraw2 <- paste(regionChrom2, altchromstart, altchromend, sep = ":" )
    }

    upper <- straw_R(sprintf("%s %s %s %s %s %i", norm, hic, regionStraw, regionStraw2, resscale, resolution))
  }

  # Full format: get symmetric data, complete missing values, replace NA's with 0's
  if(format == "full"){
    lower <- upper[ ,c(2,1,3)]
    colnames(lower) <- c("x", "y", "counts")
    combined <- unique(rbind(upper, lower))
    combinedComplete <- tidyr::complete(combined, x, y)
    combinedComplete$counts[is.na(combinedComplete$counts)] <- 0

    # Scale lower and upper bounds using zrange
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
  }

  # Sparse format: upper sparse triangular format, scale with new zrange if given
  if (format == "sparse"){
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
