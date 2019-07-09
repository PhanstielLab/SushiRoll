#' @export

gridplotHic <- function(hicdata, chrom, chromstart, chromend, palette, zrange, resolution, height = 3.25, width=3.25, x= 4.25, y= 5.5,...){

  #assuming input data is of format x, y, counts where x and y are in basepairs

  grid.newpage()
  drawpoly <- function(df, res, chrstart, chrend){
    col = rgb(df[4],df[5],df[6],maxColorValue=255)
    x = df[1]
    y = df[2]
    
    #normalized to 0 to 1 based on chromstart and chromend
    xleft = x-.5*res
    xleft.normalized = (xleft-chrstart)/(chrend-chrstart)
    xright = x+.5*res
    xright.normalized = (xright-chrstart)/(chrend-chrstart)
    ytop = y+.5*res
    ytop.normalized = (ytop - chrstart)/(chrend-chrstart)
    ybottom = y-.5*res
    ybottom.normalized = (ybottom-chrstart)/(chrend-chrstart)
    #plot all the squares
    grid.polygon(x=c(xleft.normalized,xleft.normalized,xright.normalized,xright.normalized),y=c(ybottom.normalized,ytop.normalized,ytop.normalized,ybottom.normalized),gp=gpar(col=col,fill=col))
    
  }
  
  ###SUBSET DATA###
  hicregion = hicdata[which(hicdata[,1] >= chromstart & hicdata[,1] <= chromend &
                                     hicdata[,2] >= chromstart & hicdata[,2] <= chromend),]
 
  ###CONVERT NUMBERS TO COLORS###
  
  #convert data based on zrange
  min_z = zrange[1]
  max_z = zrange[2]
  hicregion[,3][which(hicregion[,3] < min_z)] = min_z
  hicregion[,3][which(hicregion[,3] > max_z)] = max_z
  
  #use colour_values function from "colourvalues" package to convert numbers to colors
  color_vector <- colour_values(hicregion[,3], palette = palette)
  
  #append colors to hicdata and convert to rgb
  hicregion= cbind(hicregion,t(col2rgb(color_vector)))
  
  #DRAW SQUARES###
  
  # make symmetric
  hicregionrev = hicregion[,c(2,1,3,4,5,6)]
  hicregion = rbind(hicregion,setNames(hicregionrev,names(hicregion))) 
  hicregion = hicregion[!duplicated(hicregion),]
  
  
  
  if(is.null(current.vpPath()) == FALSE){
    upViewport()
  }
  #make viewport
  converted_coords = convert_coordinates(height, width, x, y)
  vp <- viewport(height=unit(height,"in"),width=unit(width,"in"),x=unit(converted_coords[1],"in"),y=unit(converted_coords[2],"in"))
  pushViewport(vp)
 
  apply(hicregion,1,drawpoly,res=resolution,chrstart = chromstart, chrend=chromend)
  
  
  #return(list(c(min_z,max_z),palette))
}

# gridplotHic(hicdata = hic.filt.full,
#          chrom = "chr8",
#          chromstart = chromstart,
#          chromend = chromend,
#          palette  = 'viridis',
#          resolution = 5000,
#          zrange=c(0,12),
#          x= 0.75,
#          y= 0.75
# )

