# RapidRicker

R Package to run spawner-recruit data quality checks and test the sensitivity of standard biological benchmarks.

Development Team: Gottfried Pestal, Charmaine Carr-Harris, Steven Cox-Rogers

**Important Notes:** 


* **WARNING**: This package is currently being fundamentally revised. For the next little while, analyses may not run,
or worse, they may run with nonsense output! Stay tuned for the next working version!** 

* This package is under development. Functions may change rapidly and substantially. Do not use these if you are not part of the development team!
* WARNING: The age comp diagnostics in the data check were not behaving properly. See [this issue](https://github.com/SOLV-Code/RapidRicker/issues/61) for details, and [this wiki page] for context. Therefore, the
corresponding steps have been commented out in the *checkSRData() fn, and the associated trigger values have been removed from the *flags_default* object.
* SR model fits and biological benchmarks are simple deterministic fits based on ```lm(logRpS ~ S)```. These are intended solely as a rapid check for the sensitivity to different ways of subsetting the data.  We envision this as a pre-screening step to plan out more formal model fitting and sensitivity testing (e.g. Bayesian models with time-varying productivity, hierarchical models with shared productivity estimates across stocks).

For now *RapidRicker* includes 3 functions that work on data for a single stock:

*  *checkSRData()*: calculates a set of diagnostics for the spawner-recruit data. Some apply to the whole series (e.g. the contrast in spawner estimates), while others flag individual observations (e.g. R/S above user-specified plausible upper bound, pointing to a potential data error in either R or S).
* *calcDetRickerBM()*: fits a simple linear regression to ```log(R/S) ~ S``` and calculates standard biological
benchmarks (Smsy, Smax, Seq, Umsy)
* *testDetRickerBM()*: runs jackknife (drop 1), retrospective (gradually add recent data) or reverse retrospective (gradually drop earlier data) variations of the simple Ricker fit.

There is also a wrapper function *RapidRicker()* that applies all of these to a data set with multiple stocks and generates a compact output object with list elements for the data check (summary tables, details), BM estimates, and sensitivity tests (e.g. BM values for each step in the retrospective, and a summary of min/max % diff of the retrospective values compared to the base case with all values).


*Update*: Basic Bayesian estimates are being added in the *calcMCMCRickerBM()* function, with the planned option of running them through the R2jags package, the Rstanarm package (with a syntax like R's lm function) or the Rstan package ("traditional" STAN). Depending on speed, these may be incorporated into the overall *RapidRicker()* function call, or stay as a standalone extension.



To get up and running with some examples, follow the *Quick Start* steps below.



## Quick Start

### Install

To install this package directly from github, use

```
install.packages("devtools") # Install the devtools package
library(devtools) # Load the devtools package.
install_github("SOLV-Code/RapidRicker", dependencies = TRUE,
                build_vignettes = FALSE)

```

**3 Warnings** will show, letting you know that the packages *R2jags*, *coda*, *rstan*, and *rstanarm* have functions with the same name: *traceplot()* and *loo()*. After installing *RapidRicker*, you need to make sure you call these functions explicitly (e.g. ```coda::traceplot()``` rather than just ```traceplot()```).


```

library(RapidRicker)	
library(tidyverse)	 
# It still needs this. 
# All the other dependencies are loaded fine...
# This seems to be by design as per
# https://www.tidyverse.org/blog/2018/06/tidyverse-not-for-packages/
# 



# check the built in data set

?SR_Sample # opens help file
head(SR_Sample) # shows the first few rows

# check the function help files

?checkSRData
?calcDetRickerBM
?testDetRickerBM

		
```




### Worked Examples


```
library(RapidRicker)
library(tidyverse)

# ---------------------------------------------
# run the wrapper function (all tests, all stocks)



# look at the default criteria for the data check 
flags_default

# run the wrapper function
rapid.ricker.out <- RapidRicker(sr_obj_m = SR_Sample, min.obs = 10,  trace=TRUE)

# check the components of the output
names(rapid.ricker.out)

# look at the data check outputs
names(rapid.ricker.out$DataCheck)
rapid.ricker.out$DataCheck$TabSeriesVal
rapid.ricker.out$DataCheck$TabSeriesFlags
rapid.ricker.out$DataCheck$TabObsFLags
head(rapid.ricker.out$DataCheck$Summary)
head(rapid.ricker.out$DataCheck$Data)

# look at the BM outputs
names(rapid.ricker.out$BM)
head(rapid.ricker.out$BM$Retro)

# look at the PercDiff outputs (sensitivity test vs. base case)
head(rapid.ricker.out$PercDiff$RetroPercDiffMin)
head(rapid.ricker.out$PercDiff$RetroPercDiffMax)


#--------------------------------------------------
# Illustration of Summary Plots

pdf(file = "RapidRicker_SamplePlots_Ranks&Tornado.pdf", onefile= TRUE, height = 8.5, width =11)

spn.df <- SR_Sample %>% select(Stock,Year,Spn) %>% 
              pivot_wider(id_cols = Year,names_from = Stock, values_from = Spn) %>% 
              select(-Year)

logRpS.df <- SR_Sample %>% select(Stock,Year,logRpS) %>% 
                  pivot_wider(id_cols = Year,names_from = Stock, values_from = logRpS) %>% 
                  select(-Year)



par(mfrow = c(1,2),mai = c(1,2,1,0.2))

# Rank stocks by Spn
plotRanking(data.df = spn.df, trim = 0,  #show full range
                         maxvars = 25, xlim = NULL, flag = "Stock5", mean.pt=FALSE)
title(main = "Ranked by Median Spn")

plotRanking(data.df = log(spn.df), trim = 0,  #show full range
            maxvars = 25, xlim = NULL, flag = "Stock5", mean.pt=FALSE)
title(main = "Ranked by Median Log(Spn)")
  
  
# Plot sensitivity range from retrospective test

smsy.range <- data.frame(Label = rapid.ricker.out$PercDiff$RetroPercDiffMin$Stock ,                      
                        Lower = rapid.ricker.out$PercDiff$RetroPercDiffMin$Smsy_p,
                        Mid = NA,
                        Upper = rapid.ricker.out$PercDiff$RetroPercDiffMax$Smsy_p)
  
plotTornado(data.df = smsy.range, keep.rank= FALSE,add.labels = TRUE,xlim = c(-100,150),
            solid.refline = 0, dashed.refline = c(-15,15))
  


umsy.range <- data.frame(Label = rapid.ricker.out$PercDiff$RetroPercDiffMin$Stock ,                      
                         Lower = rapid.ricker.out$PercDiff$RetroPercDiffMin$Umsy_p,
                         Mid = NA,
                         Upper = rapid.ricker.out$PercDiff$RetroPercDiffMax$Umsy_p)
title(main="Retrospective Test of Smsy")

plotTornado(data.df = umsy.range, keep.rank= FALSE,add.labels = TRUE,xlim = c(-100,150),
            solid.refline = 0, dashed.refline = c(-15,15))
title(main="Retrospective Test of Umsy")


dev.off()


# -----------------------------------------
# Use individual functions


data.chk <- checkSRData(SR_Sample[SR_Sample$Stock == "Stock1",])
names(data.chk)
print(data.chk$Summary)
print(head(data.chk$Data))


#  single BM calc
ricker.bm <- calcDetRickerBM(SR_Sample[SR_Sample$Stock == "Stock1",],
                              min.obs = 10)
print(ricker.bm)

# do a retrospective test (start with the first min.obs records, then add 1 at a time
ricker.retro <- testDetRickerBM(SR_Sample[SR_Sample$Stock == "Stock1",],
                               min.obs= 10,  type="retro")

print(head(ricker.retro))

# do jackknife (drop 1 record and recalculate)
ricker.jack <- testDetRickerBM(SR_Sample[SR_Sample$Stock == "Stock1",],
                                min.obs= 10,  type="jack")

print(head(ricker.jack))


# do a reverse retrospective (start with all data, and drop beginning of the series  until only min.obs are left)
ricker.revretro <- testDetRickerBM(SR_Sample[SR_Sample$Stock == "Stock1",],
                                min.obs= 10,  type="revretro")

print(head(ricker.revretro))



# Plot examples

# add a new col, for easier plotting
SR_Sample[["RpS"]] <- exp(SR_Sample$logRpS)


# single plot
spark.line(SR_Sample[SR_Sample$Stock == "Stock1",c("Year","Spn")],avg=4,type = "o",cex=0.5)
title(main = "Stock1")


# sparkline handout
pdf(file = "RapidRicker_SamplePlots_sparklines.pdf", onefile= TRUE, height = 8.5, width =11)

xlim <- range(SR_Sample$Year)


layout(matrix(c(rep(0,5),1:25),ncol=5,byrow=TRUE),heights = c(0.5,2,2,2,2,2))
par(mai = c(0.3,0.2,0.2,0.2))
#layout.show(25)

for(stk in sort(unique(SR_Sample$Stock))){
  spark.line(SR_Sample[SR_Sample$Stock == stk,c("Year","Spn")],avg=4,type = "o",cex=0.5,x.lim = xlim)
  title(main = stk)
}

title(main = "Spn", outer = TRUE, line = -2,cex.main =2, col.main = "darkblue")



layout(matrix(c(rep(0,5),1:25),ncol=5,byrow=TRUE),heights = c(0.5,2,2,2,2,2))
par(mai = c(0.3,0.2,0.2,0.2))
#layout.show(25)
rps.range <- range(SR_Sample$RpS,na.rm=TRUE)

for(stk in sort(unique(SR_Sample$Stock))){
  spark.line(SR_Sample[SR_Sample$Stock == stk,c("Year","RpS")],avg=4,type = "o",cex=0.5,x.lim = xlim,
             y.lim = c(0,25),ref.line = (c(1,5,10)))
  title(main = stk)
}

title(main = "Rec/Spn", outer = TRUE, line = -2,cex.main =2, col.main = "darkblue")


dev.off()



```



