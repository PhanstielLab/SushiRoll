convert_coordinates <- function(height, width, x, y, pageheight=11){
  
  #standardized to 8 1/2 by 11 inch piece of paper
  #given x and y coordinates in the top left
  #y coordinate based as if going from top of page
  
  ybottom = pageheight - y
  
  x1 = x + (0.5*width)
  y1 = ybottom - (0.5*height)
  
  return(list(x1,y1))
  
}
