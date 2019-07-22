#' plots HiC interaction plot
#'
#' @param hic path to .hic file or 3 column dataframe of counts
#' @param chrom chromosome of region to be plotted
#' @param chromstart chromosome start of region to be plotted
#' @param chromend chromosome end of region to be plotted
#' @param color palette to use for representing interaction scores
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param resolution the width in bp of each pixel
#' @param norm hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#' @param plottype options are "square" or "triangle"
#' @param half option for when plottype is "triangle"; options are "top" or "bottom"
#' @param raster allows for rasterization of plot, which results in quicker plotting
#' @param addlegend add legend representing scores for color range
#' @param legendlocation if addlegend == TRUE, where relative to plot to place legend; options are "right" and "top"
#' @param height height of plot
#' @param width width of plot
#' @param x x-coordinate of where to place plot
#' @param y y-coordinate of where to place plot
#'
#'
#' @export

gridplotHic <- function(hic, chrom, chromstart, chromend, palette = 'reds', zrange = NULL, resolution, norm = "NONE", plottype = "square", half = NULL, raster = TRUE, addlegend = FALSE,
                        legendlocation = NULL, height = 3.25, width=3.25, x= 0.75, y= 0.75,...){


  ###########################################################################################################################################################################################################
  ### DEFINE A FUNCTION TO PLOT SQUARES ###
  drawpoly <- function(df, resolution, chromstart, chromend){

    ## Define the color
    col = rgb(df[4], df[5], df[6], maxColorValue = 255)
    x = df[1]
    y = df[2]

    ## Get coordinates for the points of the square, normalized to 0 to 1 based on chromstart and chromend
    xleft = x - .5 * resolution
    xleft.normalized = normalize(xleft, chromstart, chromend)
    xright = x+.5*resolution
    xright.normalized = normalize(xright, chromstart, chromend)
    ytop = y+.5*resolution
    ytop.normalized = normalize(ytop, chromstart, chromend)
    ybottom = y-.5*resolution
    ybottom.normalized = normalize(ybottom, chromstart, chromend)

    ## Plot all the squares
    grid.polygon(x = c(xleft.normalized, xleft.normalized, xright.normalized, xright.normalized),
                 y = c(ybottom.normalized, ytop.normalized, ytop.normalized, ybottom.normalized), gp = gpar(col = col, fill = col))
  }

  ###########################################################################################################################################################################################################

  ### CHECK WHAT KIND OF DATA WE HAVE AND GET SUBSETTED DATAFRAME ###
  if (class(hic) == "data.frame"){
    ## Make sure dataframe is correctly subsetted and has correct header names
    if(ncol(hic) > 3){
      stop("Incorrect dataframe format.  Input a dataframe with 3 columns: x, y, counts.")
    } else {
      hicregion <- hic[which(hic[,1] >= chromstart & hic[,1] <= chromend &
                               hic[,2] >= chromstart & hic[,2] <= chromend),]
      colnames(hicregion) <- c("x", "y", "counts")

    }

    ## Scale lower and upper bounds using zrange
    if(is.null(zrange)){
      zrange <- c(0, max(hicregion$counts))
      hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
      hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
    } else {
      stopifnot(is.vector(zrange), length(zrange) == 2, zrange[2] > zrange[1])
      hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
      hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
    }

  } else if (class(hic) == "character"){
    ## Make sure the inputted file is a .hic file
    if(file_ext(hic) == "hic"){
      if(plottype == "square"){
        hicregion <- extractHiC(hic = hic, format = "full", chrom = chrom, chromstart = chromstart, chromend = chromend, resolution = resolution, zrange = zrange, norm = norm)
      } else if (plottype == "triangle"){
        hicregion <- extractHiC(hic = hic, format = "sparse", chrom = chrom, chromstart = chromstart, chromend = chromend, resolution = resolution, zrange = zrange, norm = norm)
      }
    } else {
      stop("Incorrect input. Input a .hic file or a dataframe.")
    }
  } else {

    stop("Incorrect input. Input a .hic file or a dataframe.")
  }


  ### VIEWPORTS ###

  if(is.null(current.vpPath()) == FALSE){
    upViewport()
  }
  ## Make viewport
  converted_coords = convert_coordinates(height, width, x, y)
  vp <- viewport(height=unit(height,"in"),width=unit(width,"in"),x=unit(converted_coords[1],"in"),y=unit(converted_coords[2],"in"))
  pushViewport(vp)


  ## CONVERT NUMBERS TO COLORS ##
  ## Use colour_values function from "colourvalues" package to convert numbers to colors
  color_vector <- colour_values(hicregion[,3], palette = palette)
  # Sorted color vector for use in legend
  sorted_color_vector <- colour_values(sort(hicregion[,3]), palette = palette)
  sorted_colors <- unique(sorted_color_vector)


  ## RASTERIZED PLOT ##
  if(raster == TRUE){

    ## Add color vector to hicregion dataframe
    hicregion <- cbind(hicregion,color_vector)

    ## Get the lowest color in the range to fill in matrix
    lowest_colors <- subset(hicregion, counts == min(counts))
    lowest_color <- as.character(lowest_colors[1,"color_vector"])

    ## Remove unnecessary "counts" column
    hicregion= hicregion[,c(1,2,4)]

    ## Cast dataframe into a matrix
    reshapen <- as.matrix(reshape::cast(hicregion,formula=x~y,value="color_vector"))

    if(plottype == "triangle"){

      ## Fill in all NA's with the lowest color in the range (essentially a 0)
      reshapen[is.na(reshapen)] <- lowest_color

      ## Replace the lower triangular part of the matrix with NA's to not have anything plot there
      reshapen[lower.tri(reshapen)] <- NA

      ## Reverse orientation of matrix based on columns for bottom triangle and rows for top triangle
      if (half == "bottom"){
        reshapen <- apply(reshapen,2,rev)
      } else if (half == "top"){
        reshapen <- apply(reshapen,1,rev)
      } else {
        stop("argument \"half\" is required for argument \"plottype = triangle\". ")
      }
    }


    ## Matrix already complete from extraction, reverse orientation based on columns
    if(plottype == "square"){
      reshapen <- apply(reshapen,2,rev)
    }

    grid.raster(reshapen)

  }


  ### NON-RASTERIZED PLOT ###
  if (raster == "FALSE") {

    ## Append colors to hicdata and convert to rgb
    hicregion <- cbind(hicregion,t(col2rgb(color_vector)))

    if(plottype == "triangle"){

      ## Need to get lower part of dataframe to plot bottom triangle
      if(half == "bottom"){

        lower <- hicregion[, c(2,1,3,4,5,6)]
        colnames(lower) <- c("x", "y", "counts", "red", "green", "blue")
        hicregion <- lower
      }
    }


    ## Plot squares with drawpoly function
    invisible(apply(hicregion, 1, drawpoly, resolution = resolution, chromstart = chromstart, chromend = chromend))

  }

  ### LEGEND ###

  if (addlegend == TRUE){



    min_z <- zrange[1]
    max_z <- zrange[2]

    if(legendlocation == "right"){
      legend_width <- 1/32 * width
      legend_height <- 1/4 * height
      legend_xcoord <- x + width + 0.25
      legend_ycoord <- y

      gridaddLegend(color_vector = sorted_colors, min_label = min_z, max_label = max_z, height= legend_height, width = legend_width, x = legend_xcoord, y = legend_ycoord, location = "right")

    }

    if(legendlocation == "top"){
      legend_width <- 1/4 * height
      legend_height <- 1/32 * width
      legend_xcoord <- x + (width - legend_width)/2
      legend_ycoord <- y - 0.25 - legend_height

      gridaddLegend(color_vector = sorted_colors, min_label = min_z, max_label = max_z, height= legend_height, width = legend_width, x = legend_xcoord, y = legend_ycoord, location = "top")
    }

  }


  #return sorted_color vector if choosing to add legend later
  return(sorted_colors)
}



