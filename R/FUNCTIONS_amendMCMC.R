#' amendMCMC
#' 
#'  Wrapper function that applies the calcRickerProxy() function to each MCMC sample in the calcMCMCRickerBM() output,
#'  calculates corresponding Sgen and fitted SR curves, and amends the output with summaries.
#' @param mcmc.obj output from a call to calcMCMCRickerBM() with output = "all"
#' @param sr.scale scalar used in the SR model fit (if model was fitted to Mill fish, then sr.scale = 10^6)
#' @keywords Sgen, ricker fit
#' @export



amendMCMC <- function(mcmc.obj,sr.scale =10^6){
############################

# set up the output object
mcmc.obj.out <- mcmc.obj  
  
# Extract the parameters
betas <- mcmc.obj[[1]]$MCMC$MCMC.samples[,"beta"]
kalman.check <- sum(grepl("ln.alpha\\[", dimnames(mcmc.obj[[1]]$MCMC$MCMC.samples)[[2]])) > 0
if(kalman.check){ alphas.idx <-  grepl("ln.alpha\\[", dimnames(mcmc.obj[[1]]$MCMC$MCMC.samples)[[2]]) }
if(!kalman.check){ alphas.idx <-  match("ln.alpha", dimnames(mcmc.obj[[1]]$MCMC$MCMC.samples)[[2]]) }
alphas <- mcmc.obj[[1]]$MCMC$MCMC.samples[,alphas.idx,drop=FALSE]
#alphas[alphas<0] <- NA # NEED TO DISCUSS -> now handling inside of calcRickerProxy()
num.alphas <- dim(alphas)[2]
num.mcmc <- dim(alphas)[1]

# percentiles to use for summaries
probs.use <- seq(5,95,by=5) /100

#--------------------------------------------
# SGEN CALCS

print("starting Sgen calcs")




# storage objects for sgen
sgen.quants <- matrix(NA,ncol = num.alphas, nrow=length(probs.use),dimnames = list(paste0("p",probs.use*100),
                                                                                   gsub("ln.alpha","Sgen",dimnames(alphas)[[2]])))
sgen.sample <- matrix(NA,nrow=num.mcmc,ncol = num.alphas, dimnames =list(1:num.mcmc,
                                                    gsub("ln.alpha","Sgen",dimnames(alphas)[[2]])))

# Calculate Sgen for each sample of MCMC pars (1 set for Basic Ricker and AR1 Ricker, 1 per brood year for Kalman Filter Ricker)
for(i in 1:num.alphas){	
  
  print(paste("alpha index:",i))

  vals.tmp <- mapply(calcRickerProxy, alphas[,i], betas) * sr.scale  
  check.df <- data.frame(alphas[,i], betas,vals.tmp)
  quants.tmp <- quantile(vals.tmp,probs.use,na.rm=TRUE)
  sgen.quants[,i] <- quants.tmp
  sgen.sample[,i] <- vals.tmp
  
}

# append sgen to $Percentiles Object
mcmc.obj.out[[1]]$Percentiles <- mcmc.obj.out[[1]]$Percentiles  %>%
                                  bind_rows(sgen.quants %>% t() %>% as.data.frame() %>%
                                                rownames_to_column("Variable") )

# append Sgen to $Medians object
det.a <- mcmc.obj[[1]]$Medians[mcmc.obj[[1]]$Medians$VarType == "ln_a","Det"]
det.b <- mcmc.obj[[1]]$Medians[mcmc.obj[[1]]$Medians$VarType == "b","Det"]
det.sgen <- mapply(calcRickerProxy, det.a, rep(det.b,length(det.a)))  

mcmc.obj.out[[1]]$Medians <- mcmc.obj.out[[1]]$Medians %>%
  bind_rows(
    bind_cols(
      VarType = "Sgen",
      mcmc.obj.out[[1]]$Medians %>% dplyr::filter(VarType == "Seq.c") %>% select(YrIdx,Yr),
      sgen.quants  %>% t() %>% as.data.frame() %>%
          rownames_to_column("Variable") %>% select(Variable,p10,p25,p50,p75, p90),
      Det = det.sgen) %>% 
  select (VarType, Variable, everything()) %>%
  mutate(Diff = p50 - Det ) %>%
  mutate(PercDiff = round(Diff/Det*100,1) )
  )


# append Sgen to $MCMC$MCMC.samples object
mcmc.obj.out[[1]]$MCMC$MCMC.samples <- cbind(mcmc.obj.out[[1]]$MCMC$MCMC.samples,sgen.sample)



#------------------------------------------------------------------------------
# add ricker curve percentiles to output

print("Starting fitted curve calcs")

spn.vals.use <- calcRickerProxy(a=det.a[1], b =det.b*sr.scale,  # b/c det fit currently done in fish mot scaled values 
                                spn.vals = NULL, sr.scale = sr.scale,out.type = "curve")[["spn"]]

rec.quants <- array(data = NA,dim = c(length(spn.vals.use),length(probs.use),dim(alphas)[[2]]),
                        dimnames = list(paste0("Spn",1:length(spn.vals.use)),
                                        paste0("p",probs.use*100),
                                        gsub("ln.alpha","Index",dimnames(alphas)[[2]])))

for(i in 1:num.alphas){	
print(paste("alpha index:",i))
  
exp.rec.tmp <- mapply(calcRickerProxy, a = alphas[,i], b =betas,
          MoreArgs = list(spn.vals = spn.vals.use, sr.scale = sr.scale,out.type = "rec")) %>% round()
rec.quants[,,i] <-  apply(exp.rec.tmp,MARGIN = 1,quantile, probs=probs.use,na.rm=TRUE ) %>% t()
}

mcmc.obj.out[[1]]$RickerCurve <- list(spn.vec = spn.vals.use, rec.arr.quants = rec.quants)

return(mcmc.obj.out)

} # end fn amendMCMC



