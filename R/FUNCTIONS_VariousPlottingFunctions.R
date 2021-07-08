#' Plotting Functions
#'
#' These functions create various plot components.
#' @keywords plots
#' @export


box_add <- function(perc.vec,at,width,col,bg,label){
	segments(at,perc.vec[1],at, perc.vec[5],col=col)
	rect(at - width/2, perc.vec[2], at + width/2,perc.vec[4],border=col,col = bg)
	segments(at - width/2 ,perc.vec[3],at + width/2, perc.vec[3],lwd=3,col=col,lend = 3)
	text(at,perc.vec[5],label = label,col =col,cex=0.8, adj = c(0.5,-0.3))
	}
