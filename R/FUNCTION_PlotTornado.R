#' plotTornado
#'
#' This function generates a plot showing a range for each row. Typically each row is a stock, the first col is a stock label, and the remaining columns have values for lower, mid, and upper (e.g. from a sensitivity test)
#' @param data.df  a data frame with columns Label, Lower, Mid, Upper
#' @param rank default rank is based on range between Upper and Lower. if keep.rank = TRUE, then keep the input order as is. 
#' @param add.labels if = TRUE, data.df needs another column called AddLabel which is added to the plot a labels inside the label axis (e.g. Avg abd as context)
#' @param add.pt if = TRUE, data.df needs another column called AddPt with values for an additional point to be plotted
#' @keywords sensitivity test, plot
#' @export


plotTornado <- function(data.df, add.pt = FALSE, keep.rank= FALSE, add.labels = TRUE,xlim = NULL,solid.refline = NULL,dashed.refline = NULL){


  if(!keep.rank){data.df <- data.df %>% mutate(Range = Upper - Lower) %>% arrange(-Range)	}
	
  num.units <- dim(data.df)[1]
  tick.loc <- num.units: 1#(num.units - length(data.df) + 1)

  if(is.null(xlim)){ xlim <- range(data.df[,c("Lower","Mid","Upper")],na.rm=TRUE)}


  plot(1:10, 1:10, xlim = xlim, ylim = c(1, num.units), cex.main = 0.8,
       type = "n", bty = "n", axes = FALSE, xlab = "",
       ylab = "", cex.lab = 1.5)
  
  if(!is.null(solid.refline)){abline(v=solid.refline,col="red",lty=1) }
  if(!is.null(dashed.refline)){abline(v=dashed.refline,col="red",lty=2)}

  axis(2, at = tick.loc, labels = data.df$Label,
       las = 1, lwd = 1, line = 1,cex.axis=1.1)

  if(add.labels & "AddLabel" %in% names(data.df)){
    text(par("usr")[1],tick.loc, labels = data.df$AddLabel,cex = 0.9,col="red",xpd=NA,adj=1)
    }
  x.ticks <- pretty(xlim, n = 5)

  axis(3, at = seq(-100,100,by=50))

  segments(data.df$Lower, tick.loc, data.df$Upper, tick.loc, col = "darkblue")
  points(data.df$Upper, tick.loc, pch=4, col = "darkblue",cex=0.8)
  points(data.df$Lower, tick.loc, pch=4, col = "darkblue",cex=0.8)
  points(data.df$Mid, tick.loc, pch=19, col = "darkblue")
 if(add.pt){ points(data.df$AddPt, tick.loc, pch=18, col = "darkblue")}



}