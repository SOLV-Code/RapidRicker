#' doRJAGS
#'
#' This function is a wrapper for calls to the jags() function of the R2jags package. It is designed for use inside the calcMCMCRickerBM() function, but can be used on its own if you are familiar with BUGS models. Check the code for calcMCMCRickerBM() to see the details. This function is adapted from forecasting code originally developed in collaboration with Sue Grant and Bronwyn MacDonald. Note: This requires installing JAGS from \href{https://sourceforge.net/projects/mcmc-jags/files/latest/download}{here}.
#' @param data.obj a list object with Spn and Rec (Data for 1 Stock!), as well as the required priors (details depend on model). Other variables can be there but are not used (RpS, Qual, ExpF etc)
#' @param model.fn a function that defines a BUGS model
#' @param settings a list with n.chains (2), n.burnin (20000), n.thin (60), and n.samples (50000). Default values in brackets.
#' @param inits a list of lists with inits for each chain. Depens on BUGS model
#' @param pars.track vector of text strings listing parameters to track. Depends on BUGS model
#' @param output one of "short" (only return summary stats for tracked parameters in a list object), "post" (also save posterior distribution samples to folder), or "all" (also produce pdf files with standard diagnostic plots)
#' @param out.path text string specifying  folder. if output is "post" or "all", the generated files will be stored to this folder
#' @param out.label label use in the output files if output is "post" or "all"
#' @param mcmc.seed either "default" or an integer giving the random seed for starting MCMC (R2Jags default is 123)
#' @param tracing if TRUE, diagnostic details for intermediate objects will be printed to the screen for
#' @param optional vector to provide brood years matching the data points (for labelling output)
#' @keywords R2Jags, MCMC, posterior
#' @export

doRJAGS <- function(data.obj, model.fn,
                    settings = list(n.chains=2, n.burnin=20000, n.thin=60,n.samples=50000) ,
                    inits, pars.track,
					output = "short",
					out.path = "MCMC_Out",
					out.label = "MCMC",
					mcmc.seed = "default",
					tracing = FALSE,
					sr.yrs = NA
					){





perc.vec <- seq(5,95,by=5) # %iles used in the posterior summaries


if(output == "short"){
		write.CODA <- FALSE  # write CODA txt files
		MCMC.plots <- FALSE # create MCMC diagnostic plots (traceplots etc)
		CODA.plots <- FALSE   # create plots of the posterior disctributions
		}

if(output == "post"){
		write.CODA <- FALSE  # write CODA txt files
		MCMC.plots <- FALSE # create MCMC diagnostic plots (traceplots etc)
		CODA.plots <- TRUE   # create plots of the posterior disctributions
		dir.create(out.path,showWarnings=FALSE) # creates directory, if it already exists it does nothing
		}

if(output == "all"){
		write.CODA <- TRUE  # write CODA txt files
		MCMC.plots <- TRUE # create MCMC diagnostic plots (traceplots etc)
		CODA.plots <- TRUE   # create plots of the posterior disctributions
		dir.create(out.path,showWarnings=FALSE) # creates directory, if it already exists it does nothing
		}


start.time <- proc.time()
print(paste("STARTING R2JAGS MCMC ESTIMATION FOR,", out.label, "-------------------------------------"))

if(mcmc.seed!="default"){seed.use <- mcmc.seed}
if(mcmc.seed=="default"){seed.use <- 123} # this is the default value in R2JAGS

# IMPORTANT: jags.seed argument in jags() does not work (it only applies to jags.parallel)
# Therefore need to set.seed first
set.seed(seed.use)


if("deviance" %in% pars.track){pars.track <- pars.track[pars.track != "deviance"]}

mcmc.obj <- jags(data=data.obj,
			inits=inits,
			parameters.to.save=pars.track,
			model.file=model.fn,
			DIC=TRUE,   # previusly set this to FALSE, because explicitly tracking "deviance" as one of the pars.track, but put back in for easier output
			n.chains=settings$n.chains,
			n.burnin=settings$n.burnin,
			n.thin=settings$n.thin,
			n.iter=settings$n.samples)

# if DIC = TRUE and pars.track includes "deviance", then get this error message
# In addition: Warning message:
# In FUN(X[[i]], ...) : Failed to set trace monitor for deviance
# Monitor already exists and cannot be duplicated

# BUT: if DIC = FALSE and pars.track includes "deviance", then you get a wrong DIC without error message



print(paste("MCMC - r2JAGS took",summary(proc.time()-start.time)["elapsed"]))


print(paste("STARTING OUTPUT SUMMARY FOR,", out.label, "-------------------------------------"))


# check and store current directory
base.dir <- getwd()
# output
MCMCsamples <- mcmc.obj$BUGSoutput$sims.matrix
MCMCsummary <- mcmc.obj$BUGSoutput$summary


if(tracing){
	print("Output Elements"); print(names(mcmc.obj))
	print("Model Fit"); print(mcmc.obj$model)
	print("r2jags BUGS Output Elements"); print(names(mcmc.obj$BUGSoutput))
	print("MCMC SubSample");print(MCMCsamples[1:20,]) # extract the first few rows of the chains for alpha
} # end if tracing



start.time <- proc.time()

# Save CODA in txt file (if turned on)
if(write.CODA){
			dir.create(paste(out.path,"/CODA",sep=""),showWarnings=FALSE) # creates directory, if it already exists it does nothing
			setwd(paste(out.path,"/CODA",sep=""))
			write.table(MCMCsamples,paste(out.label,"_pars.txt",sep=""))
			setwd(base.dir)
			}

# create or append an array with the MCMC samplestats
# NOTE: SEEMS THAT THESE STORAGE ARRAYS DON"T NEED TO BE EXPLICITLY REMOVED.
# THEY DISAPPEAR WHEN THE SUBROUTINE CALL ENDS BECAUSE THEY ARE NOT RETURNED TO THE PARENT FUNCTION
# SHOULD HOWEVER MAKE THIS MORE ROBUST
if(!exists("mcmc.samplestats")){
			tmp.stats <- as.array(as.matrix(MCMCsummary))
			mcmc.samplestats <- array(NA,dim=dim(tmp.stats),dimnames=list(dimnames(tmp.stats)[[1]],dimnames(tmp.stats)[[2]]))
			} # end if creating new array

# save stats from current MCMC run
mcmc.samplestats[,] <-  as.matrix(MCMCsummary) # NOTE: INCLUDES JAGS DEFAULT THINNING FOR NOW

if(tracing){ print("mcmc.samplestats");print(paste(out.label)); print(mcmc.samplestats[,])}

# create or append an array with the %iles for each tracked variable across chains
if(!exists("mcmc.percs")){
			vars.tmp <- dimnames(MCMCsamples)[[2]]
			mcmc.percs <- array(NA,dim=c(length(perc.vec),length(vars.tmp)),dimnames=list(paste("p",perc.vec,sep=""),vars.tmp))
			} # end if creating new array

mcmc.percs[,] <- apply(MCMCsamples,MARGIN=2,quantile,probs=perc.vec/100)


# create or append list object with thinned MCMC chains
if(!exists("mcmc.samples")){
		mcmc.samples <- array(NA,dim=dim(MCMCsamples),dimnames=list(1:dim(MCMCsamples)[[1]],dimnames(MCMCsamples)[[2]]))
				}
mcmc.samples[,] <- MCMCsamples


# create or append list object with DIC

if(!exists("mcmc.dic")){
			mcmc.dic <- array(NA,dim=c(1,3),dimnames=list("",c("mean(Dev)","pD","DIC")))
			} # end if creating new array

mcmc.dic[,] <- c(mcmc.samplestats["deviance","mean"],mcmc.obj$BUGSoutput$pD,mcmc.obj$BUGSoutput$DIC)

if(tracing){print("DIC ----");print(mcmc.dic[,])}


# create list object with summary tables

bugs.dic <-  data.frame(pD = mcmc.obj$BUGSoutput$pD, DIC = mcmc.obj$BUGSoutput$DIC)

bugs.summary <- mcmc.obj$BUGSoutput$summary %>%
                as.data.frame() %>% rownames_to_column("var") %>%
                dplyr::filter(!grepl("log.resid",var)) %>%
                mutate(cv = sd/mean)

# calc Gewewke diag and Gelman-Rubin Statistic (get 1 for each par)
g.score <- geweke.diag(mcmc.obj, frac1=0.1, frac2=0.5) # default setting: compare first 10% to last 50%

gelman.out <- gelman.diag(as.mcmc.list(mcmc.obj$BUGSoutput #,start = settings$n.burnin, end = numeric(0),thin = settings$n.thin
                                       ),multivariate = FALSE)
gelman.range <- range(gelman.out$psrf[,"Point est."])





names(bugs.summary) <- recode(names(bugs.summary),"2.5%" = "p2.5","25%" = "p25",
                              "50%" = "p50","75%" = "p75", "97.5%" = "p97.5")

fit.table <-   bind_cols(data.frame(pD = bugs.dic$pD,DIC = bugs.dic$DIC, max.Rhat = max(bugs.summary$Rhat),
                                    max.abs.G = max(abs(g.score$z)),
                                    min.gelman = gelman.range[1],
                                    max.gelman = gelman.range[2],
                                    med.sigma = bugs.summary %>% dplyr::filter(var=="sigma") %>% select(p50) %>% unlist(),
                                    med.deviance = bugs.summary %>% dplyr::filter(var=="deviance") %>% select(p50) %>% unlist()
                                    ),
                          bugs.summary %>% dplyr::filter(!grepl(".c",var)) %>%
							dplyr::filter(grepl("ln.alpha",var)) %>% select(cv,n.eff) %>%
								summarize(num = n(),max.cv = max(cv),
											min.n.eff = min(n.eff),
											med.n.eff = median(n.eff),
											max.n.eff = max(n.eff)) %>%
								rename_all(function(x){paste0("ln.alpha.",x)}),
                          bugs.summary %>% dplyr::filter(var=="beta") %>% select(cv,n.eff) %>%
								rename_all(	function(x){paste0("beta.",x)}),
                          bugs.summary %>% dplyr::filter(var=="sigma") %>% select(cv,n.eff) %>%
								rename_all(function(x){paste0("sigma.",x)})

                          )

beta.table <- bugs.summary %>% dplyr::filter(var=="beta") %>% select(mean,sd,cv,starts_with("p"),Rhat,n.eff) %>% select(p50,mean,cv,everything()) %>%
  dplyr::rename(median = p50)


ln.alpha.table <- bugs.summary %>%dplyr::filter(!grepl(".c",var)) %>%
							dplyr::filter(grepl("ln.alpha",var),var != "ln.alpha.prior") %>% select(mean,sd,cv,starts_with("p"),Rhat,n.eff) %>% select(p50,mean,cv,everything()) %>%
  dplyr::rename(median = p50)


if(dim(ln.alpha.table)[1] > 1){yr.idx <- 1:dim(ln.alpha.table)[1]; yrs <- sr.yrs[yr.idx]}
if(dim(ln.alpha.table)[1] == 1){yr.idx <- NA; yrs <- NA}

ln.alpha.table <- bind_cols(YrIdx = yr.idx, Year = yrs, ln.alpha.table)

if(dim(ln.alpha.table)[1] == 1){ ln.alpha.table.bookends <- ln.alpha.table }

if(dim(ln.alpha.table)[1] > 1){
	minmax.idx <- c(which.min(ln.alpha.table$median),which.max(ln.alpha.table$median))
	ln.alpha.table.bookends <- ln.alpha.table[minmax.idx,]
	}



rownames(fit.table) <- NULL
rownames(beta.table) <- NULL
rownames(ln.alpha.table) <- NULL

tables.list <-  list(fit = fit.table, beta = beta.table,
					ln.alpha.bookends = ln.alpha.table.bookends,
					ln.alpha = ln.alpha.table,
					mcmc.summary = bugs.summary)




if("S" %in% names(data.obj)){spn.tmp <- data.obj$S } 	 # need to check this: should this be na.omit(data.set[,"Spn"])

print(paste("Output processing took", summary(proc.time()-start.time)["elapsed"]))


# BUGS JAGS diagnostic plots

if (MCMC.plots){

start.time <- proc.time()

print(paste("STARTING BUGS/JAGS DIAGNOSTICS FOR,", out.label, "-------------------------------------"))
# NOTE this calculates some diagnostics, and creates a pdf of plots if plotting is turned on

dir.create(paste(out.path,"/PLOTS",sep=""),showWarnings=FALSE) # creates directory, if it already exists it does nothing

pdf(paste(out.path,"/PLOTS/", paste(out.label,"DiagnosticPlots_JAGS.pdf",sep="_"),sep=""),width=8.5, height=8.5, onefile=TRUE) ; par(mfrow=c(1,1))  # change dir and start pdf

plot(mcmc.obj)# basic plot

# plot.jags does not include a density plot like the densplot() in the coda package
# could just do a hist() here?

R2jags::traceplot(mcmc.obj,ask=FALSE)# traceplot() not in r2OpenBUGS

dev.off(); setwd(base.dir)  # close pdf and return to working folder

print(paste("JAGS diagnostic plots took", summary(proc.time()-start.time)["elapsed"]))

} # end if  JAGS.diag.plots



# OUTPUT - CODA

if (CODA.plots){

start.time <- proc.time()

print(paste("STARTING CODA DIAGNOSTICS FOR,", paste(out.label), "-------------------------------------"))
# NOTE this calculates some diagnostics, and creates a pdf of plots if plotting is turned on

dir.create(paste(out.path,"/PLOTS",sep=""),showWarnings=FALSE) # creates directory, if it already exists it does nothing

pdf(paste(out.path,"/PLOTS/", paste(out.label,"DiagnosticPlots_CODA.pdf",sep="_"),sep=""),width=8.5, height=8.5, onefile=TRUE) ; par(mfrow=c(1,1))  # change dir and start pdf
print("starting conversion to coda file")

# convert output to make usable for diagnostics from coda package
coda.obj1 <- as.mcmc(mcmc.obj$BUGSoutput$sims.matrix)

print("conversion to coda file successful")

#xyplot(coda.obj1)  # -> not creating any plots WHY?
plot(coda.obj1)
#gelman.plot(coda.obj2)  # NOT WORKING YET
crosscorr.plot(coda.obj1,main="crosscorr.plot")
cumuplot(coda.obj1)
densplot(coda.obj1)
#print("flag: starting geweke plot")
#warning("SKIPPING GEWEKE PLOT FOR NOW in mcmc.sub()")
#geweke.plot.MOD(coda.obj1)

dev.off(); setwd(base.dir)  # close pdf and return to working folder

print(paste("CODA diagnostic plots took", summary(proc.time()-start.time)["elapsed"]))

} # end if  CODA.diag.plots


#############################################################
print("CREATING OUTPUT OBJECT -------------------------------------")


# CREATING OUTPUT LIST OBJECT (ONLY PARTLY IMPLEMENTED FOR NOW)
out.list <- list(mcmc.call=out.label,mcmc.settings=unlist(settings))

if(output %in% c("short","post","all")){out.list<-c(out.list,list(SampleStats=mcmc.samplestats, MCMC.Percentiles=mcmc.percs,Conv.Info="TBI",
					DIC=mcmc.dic, tables = tables.list))}

if(output %in% c("post","all")){out.list<-c(out.list,list(Data=data.obj))}

if(output %in% c("post","all")){out.list<-c(out.list,list(MCMC.samples=mcmc.samples))}

if(output =="all"){out.list<-c(out.list,list(MCMC.obj=mcmc.obj))}

print(names(out.list))


return(out.list)

} #end doRJAGS











