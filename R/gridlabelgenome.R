#' @export

gridlabelgenome <- function(chrom, chromstart, chromend, scale="bp",width = 3.25, x = 0.75, y = 4, fontsize = 10){
  require(grid)
  
  if (scale == "bp"){
    fact = 1
  }
  
  if (scale == "Mb"){
    fact = 1000000
  }
  
  if (scale == "Kb"){
    fact = 1000
  }
  
  chromstartlabel = chromstart/fact
  chromendlabel = chromend/fact

  ##PLOT THE AXIS##

  if(is.null(current.vpPath()) == FALSE){
    upViewport()
  }
  
  converted_coords = convert_coordinates(height=0.125, width, x, y)
  
  vp <- viewport(width=unit(width,"in"),height = unit(.125,"in"),x=unit(converted_coords[1],"in"),y=unit(converted_coords[2],"in"))
  pushViewport(vp)
  grid.segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp=gpar(col="grey"))
  #grid.xaxis(at=c(0,1),label=c(paste(chromstartlabel,scale,sep=" "),paste(chromendlabel,scale,sep=" ")))
  grid.text(chrom, x=0.5, y=0, gp=gpar(fontface="bold",fontsize = fontsize,col="grey"))
  grid.text(paste(chromstartlabel,scale,sep=" "), x= 0, y=0, just="left",gp=gpar(fontsize = fontsize,col="grey"))
  grid.text(paste(chromendlabel,scale,se=" "), x=1.025,y=0,just="right",gp=gpar(fontsize = fontsize,col="grey"))


} 





