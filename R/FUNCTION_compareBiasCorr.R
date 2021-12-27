#' compareBiasCorr
#'
#' This function generates a summary of raw and bias-corrected  values for SR parameters and resulting biological
#' benchmarks, using 
#' @param bm bm_obj output from \code{\link[RapidRicker]{calcMCMCRickerBM}}.
#' @return a data frame with 1 column for each parameter (ln.a,Seq,Smsy,Sgen) and different versions in rows. Versions include: deterministic estimate with and without bias correction. Mean estimate of MCMC posteriors with (1) no bias correction, (2) bias correction using mean ln.a and mean sigma, and (3) bias correction applied to each MCMC sample and then calculate the mean. Median estimate and perentiles of MCMC posteriors with (1) no bias correction, (2)  bias correction applied to each MCMC sample.
#' @keywords bias correction, MCMC, posterior
#' @export


compareBiasCorr <- function(bm_obj){


means.df <- apply(bm_obj$MCMC,MARGIN = 2,mean) %>% as.data.frame()
names(means.df) <- "Mean"


bm.obj.summary <- bm.out$Summary %>% dplyr::rename(Median = p50) %>%
            left_join(means.df %>% rownames_to_column("Variable"),by="Variable")


vars.vec <- c("ln.alpha", "Seq","Smsy","Sgen","SgenRatio")

# no bias correction
out.none <- bm.obj.summary %>% dplyr::filter(VarType  %in% c("beta","sigma",vars.vec)) %>%
  select(Variable,Det,Mean,Median, p10,p25,p75,p90) %>%
  column_to_rownames("Variable") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("Stat") %>%
  mutate(BiasCorr = "None") %>%
  select(Stat, BiasCorr,everything())



# bias corr using mean ln.a;pha and mean sigma

mean.corr <- out.none %>% dplyr::filter(out.none$Stat == "Mean" & out.none$BiasCorr == "None" ) %>%
                select(Stat, BiasCorr, beta, ln.alpha, sigma)
mean.corr$BiasCorr <- "Mean"
mean.corr$ln.alpha <- mean.corr$ln.alpha + mean.corr$sigma^2/2
mean.corr

mean.corr <- calcRickerOtherBM(X = mean.corr,out.type="Full") %>% select(-Smax)
mean.corr <- calcRickerSmsy(X = mean.corr , method = bm_obj$methods$Smsy ,sr.scale = sr.scale.use, out.type = "Full")
mean.corr <- calcRickerSgen(X = mean.corr , method = bm_obj$methods$Sgen,sr.scale = sr.scale.use, out.type = "Full") %>%
              select(-SmsyCalc,-SgenCalc) %>% dplyr::rename(SgenRatio = Ratio)
mean.corr



# bias corr applied to each MCMC sample

out.cs<- bm.obj.summary %>% dplyr::filter(VarType  %in% c("beta","sigma",paste0(vars.vec,".c"))) %>%
  select(Variable,Det,Mean,Median, p10,p25,p75,p90) %>%
  column_to_rownames("Variable") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("Stat") %>%
  mutate(BiasCorr = "Sample") %>%
  select(Stat, BiasCorr,everything())

names(out.cs) <- gsub(".c","",names(out.cs))
out.cs[out.cs$Stat == "Det" & out.cs$BiasCorr == "Sample" ,"BiasCorr"] <- "Mean"



return(bind_rows(out.none,mean.corr,out.cs) %>% arrange(Stat))


} # end testBiasCorr