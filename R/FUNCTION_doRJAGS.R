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
#' @param tracing if TRUE, diagnostic details for intermediate objects will be printed to the screen for debugging
#' @keywords R2Jags, MCMC, posterior
#' @export

doRJAGS <- function(data.obj, model.fn, 
                    settings = list(n.chains=2, n.burnin=20000, n.thin=60,n.samples=50000) ,
                    inits, pars.track,
					output = "short",
					out.path = "MCMC_Out",
					out.label = "MCMC",
					mcmc.seed = "default",
					tracing = FALSE
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

mcmc.obj <- jags(data=data.obj, 
			inits=inits, 
			parameters.to.save=pars.track, 
			model.file=model.fn,  
			DIC=TRUE,   # set this to FALSE, because explixitly tracking "deviance" as one of the pars.track
			n.chains=settings$n.chains, 
			n.burnin=settings$n.burnin, 
			n.thin=settings$n.thin, 
			n.iter=settings$n.samples)
	
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
					DIC=mcmc.dic))}
				
if(output %in% c("post","all")){out.list<-c(out.list,list(Data=data.obj))}				
				
if(output %in% c("post","all")){out.list<-c(out.list,list(MCMC.samples=mcmc.samples))}

if(output =="all"){out.list<-c(out.list,list(MCMC.obj=mcmc.obj))}

print(names(out.list))


return(out.list) 

} #end doRJAGS











