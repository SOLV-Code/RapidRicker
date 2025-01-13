library(RapidRicker)
library(tidyverse)



?SR_Sample # opens help file
head(SR_Sample) # shows the first few rows

# check the function help files

?checkSRData
?calcDetRickerBM
?testDetRickerBM




#look at the default criteria for the data check
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





