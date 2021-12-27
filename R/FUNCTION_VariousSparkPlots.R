#' Sparklines
#'
#' These functions generate various minimalist plots that can be combined into a multipanel dashboard. Check the worked examples in \href{https://github.com/SOLV-Code/RapidRicker}{Package ReadMe} for illustrations.
#' @param x a data frame with Year and Some Value (for spark.line) or a data frame with Year and lower,mid,upper (for spark.band) or a text string (for spark.text). 
#' @param x.lim plot limits for x-axis (Year). Default is NULL, producing default R axis range.
#' @param y.lim plot limits for y-axis (Value).Default is NULL, producing default R axis range.
#' @param avg number of years used for rolling avg (applies to spark.line)
#' @param ref.line  1 or more values for a horizontal reference line (applies to spark.line and spark.band)
#' @param ref.line.v 1 or more values for a vertical reference line (applies to spark.line and spark.band)
#' @param ref.band   vector with 2 values giving a horizontal band
#' @param type plot type. default = "l" for lines. other options are ("p" for points, "o" for points and lines)
#' @param pch plotting character (19 = solid circle)
#' @param cex = point expansion. default = 0.2 (20 percent of default)
#' @keywords sparkline
#' @export
#' @examples
#' spark.line(SR_Sample[SR_Sample$Stock == "Stock1","Spn"])


spark.line <- function(x,avg = 4, x.lim = NULL, y.lim = NULL, ref.line = NULL, ref.line.v = NULL, ref.band = NULL, type="l", pch=19,cex=0.2){


n.vals <- sum(!is.na(x[,2]))

if(n.vals >0){ #plot only if have data


# fill in missing years in input with NA (to avoid plotting issues)

yrs.vec <- min(x$Year):max(x$Year)
missing.yrs <- yrs.vec[!(yrs.vec %in% x$Year)] 

if(length(missing.yrs)>0){

yrs.fill <- data.frame(Year = missing.yrs, Value = NA) 
names(yrs.fill) <- names(x)

x <- rbind(x,yrs.fill) %>% arrange(Year)

#print(x)

}


if(is.null(x.lim)){x.lim <- range(x[,1],na.rm=TRUE)}
if(is.null(y.lim)){
            if(sum(!is.na(x[,2]))>0){y.lim <- range(x[,2],na.rm=TRUE)}
            if(sum(!is.na(x[,2]))==0){y.lim <- c(0,5)}
            }


plot(x[,1],x[,2],axes=FALSE,bty="n",type=type,col="darkblue",pch=pch,cex=cex,
     xlim = x.lim,ylim=y.lim, xlab="",ylab="",lwd=2)


if(!is.null(ref.band)){
			rect(par("usr")[1],ref.band[1],par("usr")[2],ref.band[2],col="lightgrey",border="lightgrey")
			
			}

if(!is.null(ref.line)){abline(h=ref.line,col="darkgrey",lwd=2,lty=1)}

if(!is.null(ref.line.v)){abline(v=ref.line.v,col="darkgrey",lwd=2,lty=1)}



if(!is.null(ref.line) | !is.null(ref.line.v) ){
		lines(x[,1],x[,2],type=type,col="darkblue",pch=pch,cex=cex,lwd=2)
		}



if(!is.null(avg)){
  if(n.vals >avg){
  avg.vals <- stats::filter(x[,2],filter = rep(1/avg,avg),sides = 1)
                  lines(x[,1],avg.vals,col="red",lwd=3)
  }}

} # end if doing plot


if(n.vals == 0){ # empty plot
  plot(1:5,1:5,axes=FALSE,bty="n",type="n",xlab="",ylab="")
}


} # end spark.line




#' @rdname spark.line
#' @export

spark.band <- function(x,avg = 4,x.lim = NULL,y.lim = NULL,ref.line = NULL){

  n.vals <- sum(!is.na(x[,2]))

  if(n.vals >0){ # only do plot if have data

  if(is.null(x.lim)){x.lim <- range(x[,1],na.rm=TRUE)}
  if(is.null(y.lim)){
    if(sum(!is.na(x[,2]))>0){y.lim <- range(x[,2:4],na.rm=TRUE)}
    if(sum(!is.na(x[,2]))==0){y.lim <- c(0,5)}
  }

  plot(x[,1],x[,3],axes=FALSE,bty="n",type="l",col="darkblue",
       xlim = x.lim,ylim=y.lim, xlab="",ylab="")


  polygon(c(x[,1],rev(x[,1])),
          c(x[,4],rev(x[,2])),
          col="lightgrey",border=NA)

  if(!is.null(ref.line)){abline(h=ref.line,col="darkgrey",lwd=2,lty=2)}


  lines(x[,1],x[,3],type="l",lwd=4,col="darkblue")
  } # end if doing plot

  if(n.vals == 0){ # empty plot
    plot(1:5,1:5,axes=FALSE,bty="n",type="n",xlab="",ylab="")
  }


}


#' @rdname spark.line
#' @export

spark.text<- function(x,pos=c(5,5),adj=0.5,cex =1,font = 1,col="black"){
# x is a text string
# pos specifies text position on a 10x10 grid
# cex,adj =  as per text()
plot(1:10,1:10,bty="n",axes=FALSE,bty="n",type="n",xlab="",ylab="")
text(pos[1],pos[2],labels=x,cex=cex,adj=adj,xpd=NA,font=font,col=col)
}

#' @rdname spark.line
#' @export

spark.sep <- function(col="lightblue",lty=1,lwd=2){
# plots a user-specified horizontal line in the middle of the panel
    plot(1:10,1:10,bty="n",axes=FALSE,bty="n",type="n",xlab="",ylab="")
  abline(h=5,col=col,lty=lty,lwd=lwd)

}



