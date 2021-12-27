#' plotRanking
#'
#' This function generates a plot showing range and medians for columns in a table, ranked by median or mean.  Original function created for the \href{https://github.com/SOLV-Code/SPATFunctions-Package}{SPATFunctions Package}, but adapted here to customize formatting for report figures and avoid adding a dependency.
#' @param data.df  a data frame where each column is a set of values. Column labels become plot labels
#' @param trim an integer value between 0 and 49 to specify the percentile to be plotted for the whiskers (e.g. trim = 25 plots upper/lower quartiles). 
#' @param maxvars integer specifying the maximum number of items  to plot (e.g. 30 = plot top 30)
#' @param xlim vector to specify the range for the horizontal axis
#' @param flag  optional text string specifying a column to highlight inthe ranked plot (e.g. a Stock Name)
#' @param mean.pt if TRUE, show mean instead of median. 
#' @keywords rank, plot
#' @export


plotRanking <- function (data.df, trim = 25, 
maxvars = NULL, xlim = NULL, flag = NULL,mean.pt=FALSE){
    options(scipen = 100000)
    if (length(data.df) == 1) {
        plot(1:5, 1:5, type = "n", xlab = "", ylab = "", 
            axes = FALSE)
        text(2.5, 2.5, "Need at least 2 series for this plot")
		
		data.summary <- NA
		
    }
	
	
    if (length(data.df) > 1) {
	
	
        data.summary <- as.data.frame(t(apply(data.df, MARGIN = 2, 
            FUN = quantile, probs = seq(0, 1, by = 0.01), na.rm = TRUE)))
        
       
		data.summary <- data.summary[, c(paste0(trim, "%"), 
            "50%", paste0(100 - trim, "%"))]
        
		if(mean.pt){data.summary[,"50%"] <- colMeans(data.df,na.rm=TRUE)}
		
		data.summary <- data.summary[order(data.summary[, "50%"], 
            decreasing = TRUE), ]			

		if (is.null(xlim)) {xlim <- range(data.summary,na.rm=TRUE)}
        if (!is.null(maxvars)){num.units <- maxvars}
        if (is.null(maxvars)){num.units <- length(data.df)}
		
        tick.loc <- num.units:(num.units - length(data.df) + 1)
    

    plot(1:10, 1:10, xlim = xlim, ylim = c(1, num.units), cex.main = 0.8, 
        type = "n", bty = "n", axes = FALSE, xlab = "", 
        ylab = "", cex.lab = 1.5)
    axis(2, at = tick.loc, labels = dimnames(data.summary)[[1]], 
        las = 1, lwd = 1, line = 1,cex.axis=1.1)
    x.ticks <- pretty(xlim, n = 5)
    x.ticks <- x.ticks[x.ticks < xlim[2] & x.ticks > xlim[1]]
    if (max(x.ticks) <= 10^3) {
        x.ticks.labels <- x.ticks
    }
    if (max(x.ticks) > 10^3) {
        x.ticks.labels <- paste(x.ticks/10^3, "k", sep = "")
    }
    if (max(x.ticks) > 10^4) {
        x.ticks.labels <- paste(x.ticks/10^4, "*10k", sep = "")
    }
    if (max(x.ticks) > 10^5) {
        x.ticks.labels <- paste(x.ticks/10^5, "*100k", 
            sep = "")
    }
    if (max(x.ticks) > 10^6) {
        x.ticks.labels <- paste(x.ticks/10^6, "M", sep = "")
    }
    axis(3, at = x.ticks, labels = x.ticks.labels)
	
    if (!is.null(flag)) {
	
        flag.idx <- dimnames(data.summary)[[1]] == flag
        abline(h = c(tick.loc[flag.idx] + 0.5, tick.loc[flag.idx] - 
            0.5), col = "tomato", lty = 2, xpd = TRUE)
    }
	
    lines(data.summary[, "50%"], tick.loc, type = "b", 
        pch = 19, cex = 1.6, col = "darkblue", xpd = NA)
    segments(data.summary[, paste0(trim, "%")], tick.loc, 
        data.summary[, paste0(100 - trim, "%")], tick.loc, 
        col = "darkblue")
		
		
	} # end if doing plots	
	
return(data.summary)
	

}