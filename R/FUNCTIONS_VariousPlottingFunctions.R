#' Plotting Functions
#'
#' These functions create various plot components.
#' @keywords plots
#' @export


box_add <- function (perc.vec, at, width, col, bg, label,dir="v") {
  
  if(dir == "v"){
    
    segments(at, perc.vec[1], at, perc.vec[5], col = col,xpd=NA)
    rect(at - width/2, perc.vec[2], at + width/2, perc.vec[4], 
         border = col, col = bg,xpd=NA)
    segments(at - width/2, perc.vec[3], at + width/2, perc.vec[3], 
             lwd = 3, col = col, lend = 3,xpd=NA)
    
  }
  
  if(dir == "h"){
    
    segments(perc.vec[1], at, perc.vec[5], at, col = col,xpd=NA)
    rect(perc.vec[2], at - width/2,   perc.vec[4], at + width/2,
         border = col, col = bg,xpd=NA)
    segments( perc.vec[3], at - width/2, perc.vec[3], at + width/2,
              lwd = 3, col = col, lend = 3,xpd=NA)
    
  }  
 }
 
 
 
 


kernel_margin <- function(x,at, width, col,dir="v"){
  
dens.calc <- density(x)
  #head(a.margin)
  #plot(a.margin)
  
dens.calc$y <- at + dens.calc$y/max(dens.calc$y)*width
if(dir == "v"){lines(x = dens.calc$y, y = dens.calc$x,col=col,xpd=NA,lwd=2)}
if(dir == "h"){lines(x = dens.calc$x, y = dens.calc$y,col=col,xpd=NA,lwd=2)}
  
}
