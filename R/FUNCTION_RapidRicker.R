#' RapidRicker
#'
#' This function is a wrapper for testDetRickerBM(). It loops throug multiple stocks and applies all the sensitivity test variations, then calculates the perc differences to base case. For calculation details, refer to calcDetRickerBM and testDetRickerBM.
#' @param sr_obj_m a data frame with Stock, Year and Spn, logRpS (Data for multiple stocks!). Other variables can be there but are not used (RpS, Qual, ExpF etc). 
#' @param min.obs min number of S-R pairs needed to fit a model
#' @param trace if TRUE, print various intermediate diagnostic output to the console
#' @param flags a data frame with Label, Lower, Upper. Lower and Upper define the triggers for values to be flagged. Labels have to match the criteria used in the function(LIST). If flags = NULL, then it uses the built-in object flags_default.
#' @param data.check if TRUE, do the data check (either data.check or bm.test or both must be true!)
#' @param bm.test  if TRUE, do the benchmark sensitivity tests (either data.check or bm.test or both must be true!)
#' @keywords sensitivity test
#' @export
#' @examples
#' rapid.ricker.out <- RapidRicker(SR_Sample,min.obs = 10,  trace=FALSE)


RapidRicker <- function(sr_obj_m,min.obs = 10, trace = TRUE, flags = NULL, data.check = TRUE, bm.test = TRUE,sr.scale=10^6){

# NEED TO FIX:
# - how trace = TRUE is handled throughout the subroutine calls
# - reduce what's currently printed to screen in the subroutines.

stk.list <- sort(unique(sr_obj_m$Stock))

  print(paste(length(stk.list), "Stocks"))
  print(stk.list)


# --------------------------------------------------------------------------
# Part 0: Run the data check

if(data.check){

if(is.null(flags)){flags <- flags_default}

# output objects
data.check.summary <- list(NULL)
data.check.data <- list(NULL)


for(stk in unique(sr_obj_m$Stock)){
  
  print("----------------------------")
  print(stk)
  
  # prep the data
  sr.sub <- sr_obj_m %>% dplyr::filter(Stock == stk) 
  
  if(dim(sr.sub)[1] > 0 ){
    sr.in <-  sr.sub[,names(sr.sub) %in% c("Year","Spn","Rec","SpnExp") | grepl("RecAge",names(sr.sub))]
    rec.age.idx <- grepl("RecAge",names(sr.in))
    # convert to Proportion (in case it isn't), then round
    sr.in[,rec.age.idx] <- round(sr.in[,rec.age.idx]/ rowSums(sr.in[,rec.age.idx],na.rm=TRUE),4) 
    
    data.check.tmp <- checkSRData(sr_obj = sr.in, flags  = flags)
    
    data.check.summary[[stk]] <-  cbind(Stock = stk, data.check.tmp$Summary)
    data.check.data[[stk]] <-  cbind(Stock = stk, data.check.tmp$Data)
  }
  
  
} #end looping through stocks

# convert list to data frame
data.check.summary <-  bind_rows(data.check.summary)
data.check.data <-  bind_rows(data.check.data)
  
  
table.series.val <- data.check.summary %>% dplyr::filter(Scope == "Series") %>% 
          pivot_wider(id_cols = Stock, names_from = Label, values_from = MetricVal)
            
table.series.flags <- data.check.summary %>% dplyr::filter(Scope == "Series") %>% 
  pivot_wider(id_cols = Stock, names_from = Label, values_from = Flagged)

table.obs.numflagged <- left_join(
  table.series.val %>% select(Stock, NumObs),   
  data.check.summary %>% dplyr::filter(Scope == "Obs") %>% 
  pivot_wider(id_cols = Stock, names_from = Label, values_from = NumFlagged),
  by = "Stock")
  
  
  

data.check.list <- list(TabSeriesVal = table.series.val,
						TabSeriesFlags = table.series.flags,
						TabObsFLags = table.obs.numflagged,
						Summary = data.check.summary, Data = data.check.data)


} # end if data.check


if(bm.test){

#-----------------------------------------------------------------------
# Part 1: Simple Ricker BM

if(trace){print("starting simple Ricker BM")}

# OLD
#bm.cols <- c("n_obs", "ln_a","ln_a_c","a","b","sd","Smax","Seq","Seq.c","Smsy","Smsy.c",
#             "Umsy","Umsy.c")

# NEW
bm.cols <- c("n_obs","ln.alpha","ln.alpha.c","alpha","alpha.c","beta","sigma","Seq",
"Smax","Seq.c","Smsy","Smsy.c","Umsy","Umsy.c","Sgen","Sgen.c","SgenRatio","SgenRatio.c")


bm.det.store <- as.data.frame(matrix(NA,nrow=length(stk.list),ncol= length(bm.cols),
                                     dimnames=list(stk.list,bm.cols)))


for(stk in stk.list){

  print("------------")
  print(stk)
  sr.sub <- sr_obj_m %>% dplyr::filter(Stock == stk)
  if(trace){print(sr.sub)}
  
  fit.tmp <- calcDetModelFit(sr_obj = sr.sub,
  sr.scale = sr.scale, min.obs=min.obs,resids = FALSE, fn.use = "lm", ar1 = FALSE)
   bm.tmp <- calcDetRickerBM(fit_obj = fit.tmp,sr.scale = sr.scale)
  if(trace){print(bm.tmp)}
  
  bm.det.store[stk,] <- bm.tmp
}

bm.det.store <- bm.det.store %>% rownames_to_column(var="Stock")
if(trace){print(head(bm.det.store))}


#-----------------------------------------------------------------------
# Part 2: Drop-1 Jackknife Simple Ricker BM


jack.bm.det.store <- bm.det.store[FALSE,] %>% mutate(DropYr = NA) %>% select(Stock, DropYr, everything())

for(stk in stk.list){

  print("------------")
  print(stk)

  sr.sub <- sr_obj_m %>% dplyr::filter(Stock == stk) 

  if(sum(!is.na(sr.sub$Spn) & !is.na(sr.sub$logRpS) ) > min.obs){

    jack.tmp <- testDetRickerBM(sr_obj = sr.sub,sr.scale = sr.scale,min.obs=min.obs,type="jack")
    jack.tmp <- jack.tmp %>% mutate(Stock = stk) %>% select(Stock, DropYr,everything())
    jack.bm.det.store <- rbind(jack.bm.det.store,jack.tmp)
    }

}

if(trace){print(head(jack.bm.det.store))}



#-----------------------------------------------------------------------
# Part 3: Retrospective

retro.bm.det.store <- bm.det.store[FALSE,] %>% mutate(DropYr = NA) %>% select(Stock, DropYr, everything())

for(stk in stk.list){

  print("------------")
  print(stk)

  sr.sub <- sr_obj_m %>% dplyr::filter(Stock == stk) 
  
  if(sum(!is.na(sr.sub$Spn) & !is.na(sr.sub$logRpS) ) > min.obs){

    retro.tmp <- testDetRickerBM(sr_obj = sr.sub,sr.scale = sr.scale,min.obs=min.obs,type="retro")
    retro.tmp <- retro.tmp %>% mutate(Stock = stk) %>% select(Stock, UpToYr,everything())
    retro.bm.det.store <- rbind(retro.bm.det.store,retro.tmp)
  }

}


if(trace){print(head(retro.bm.det.store))}




#-----------------------------------------------------------------------
# Part 4: Reverse Retrospective

revretro.bm.det.store <- bm.det.store[FALSE,] %>% mutate(SinceYr = NA) %>% select(Stock, SinceYr, everything())

for(stk in stk.list){

  print("------------")
  print(stk)

  sr.sub <- sr_obj_m %>% dplyr::filter(Stock == stk) 
  
  if(sum(!is.na(sr.sub$Spn) & !is.na(sr.sub$logRpS) ) > min.obs){

    revretro.tmp <- testDetRickerBM(sr_obj = sr.sub,sr.scale = sr.scale,min.obs=min.obs,type="revretro",trace=FALSE)
    revretro.tmp <- revretro.tmp %>% mutate(Stock = stk) %>% select(Stock, SinceYr,everything())
    revretro.bm.det.store <- rbind(revretro.bm.det.store,revretro.tmp)
  }

}



if(trace){print(head(revretro.bm.det.store))}

#-----------------------------------------------------------------------
# Part 5: Drop 2 largest Spn (to check how much "pull" they have on the fit)


drop2.bm.det.store <- as.data.frame(matrix(NA,nrow=length(stk.list),ncol= length(bm.cols),
                                           dimnames=list(stk.list,bm.cols)))

for(stk in stk.list){

  print("------------")
  print(stk)



  sr.sub <- sr_obj_m %>% dplyr::filter(Stock == stk) 

    # drop the 2 R,S pairs with the largest Spn
    drop.val <- sort(sr.sub$Spn[!is.na(sr.sub$logRpS)],decreasing = TRUE)[2]
    sr.sub <- sr.sub %>% dplyr::filter(Spn < drop.val)


  fit.tmp <- calcDetModelFit(sr_obj = sr.sub,
  sr.scale = sr.scale, min.obs=min.obs,resids = FALSE, fn.use = "lm", ar1 = FALSE)
  
  bm.tmp <- calcDetRickerBM(fit_obj = fit.tmp,sr.scale = sr.scale)
  if(trace){print(bm.tmp)}

  drop2.bm.det.store[stk,] <- bm.tmp

}

drop2.bm.det.store <- drop2.bm.det.store %>% rownames_to_column(var="Stock")
if(trace){print(head(drop2.bm.det.store))}



bm.list <- list(BaseCase = bm.det.store,
                Retro= retro.bm.det.store,
                RevRetro = revretro.bm.det.store,
                Jack = jack.bm.det.store,
                Drop2 = drop2.bm.det.store)



# calculate the perc differences to base case

# Jackknife
jack.max <- left_join(bm.list$BaseCase %>% select(Stock),
                      bm.list$Jack %>% select(-DropYr, -n_obs) %>% group_by(Stock) %>% summarize_all(list(max=max),na.rm=TRUE),
                      by="Stock") # need to merge it with base case to get full list of stocks, with NA for those not in jackknife output
names(jack.max) <- gsub("_max","",names(jack.max))

jack.min <- left_join(bm.list$BaseCase %>% select(Stock),
                      bm.list$Jack %>% select(-DropYr, -n_obs) %>% group_by(Stock) %>% summarize_all(list(min=min),na.rm=TRUE),
                      by="Stock")
names(jack.min) <- gsub("_min","",names(jack.min))

jack.max.effect <- cbind(Stock = bm.list$BaseCase[,1],
                         round((jack.max[,-1] - bm.list$BaseCase[,-c(1,2)])/bm.list$BaseCase[,-c(1,2)]*100))

jack.min.effect <- cbind(Stock = bm.list$BaseCase[,1],
                         round((jack.min[,-1] - bm.list$BaseCase[,-c(1,2)])/bm.list$BaseCase[,-c(1,2)]*100))


# Retro

retro.max <- left_join(bm.list$BaseCase %>% select(Stock),
                       bm.list$Retro %>% select(-UpToYr, -n_obs) %>% group_by(Stock) %>% summarize_all(list(max=max),na.rm=TRUE),
                       by="Stock") # need to merge it with base case to get full list if stocks, with NA for those not in retro output
names(retro.max) <- gsub("_max","",names(retro.max))

retro.min <- left_join(bm.list$BaseCase %>% select(Stock),
                       bm.list$Retro %>% select(-UpToYr, -n_obs) %>% group_by(Stock) %>% summarize_all(list(min=min),na.rm=TRUE),
                       by="Stock") # need to merge it with base case to get full list if stocks, with NA for those not in retro output
names(retro.min) <- gsub("_min","",names(retro.min))

retro.max.effect <- cbind(Stock = bm.list$BaseCase[,1],
                          round((retro.max[,-1] - bm.list$BaseCase[,-c(1,2)])/bm.list$BaseCase[,-c(1,2)]*100))

retro.min.effect <- cbind(Stock = bm.list$BaseCase[,1],
                          round((retro.min[,-1] - bm.list$BaseCase[,-c(1,2)])/bm.list$BaseCase[,-c(1,2)]*100))


# Reverse Retro

# calculate the % change for largest and smallest reverse retrospective outputs

revretro.max <- left_join(bm.list$BaseCase %>% select(Stock),
                          bm.list$RevRetro %>% select(-SinceYr, -n_obs) %>% group_by(Stock) %>% summarize_all(list(max=max),na.rm=TRUE),
                          by="Stock") # need to merge it with base case to get full list if stocks, with NA for those not in retro output
names(revretro.max) <- gsub("_max","",names(revretro.max))

revretro.min <- left_join(bm.list$BaseCase %>% select(Stock),
                          bm.list$RevRetro %>% select(-SinceYr, -n_obs) %>% group_by(Stock) %>% summarize_all(list(min=min),na.rm=TRUE),
                          by="Stock") # need to merge it with base case to get full list if stocks, with NA for those not in retro output
names(revretro.min) <- gsub("_min","",names(revretro.min))



revretro.max.effect <- cbind(Stock = bm.list$BaseCase[,1],
                             round((revretro.max[,-1] - bm.list$BaseCase[,-c(1,2)])/bm.list$BaseCase[,-c(1,2)]*100))

revretro.min.effect <- cbind(Stock = bm.list$BaseCase[,1],
                             round((revretro.min[,-1] - bm.list$BaseCase[,-c(1,2)])/bm.list$BaseCase[,-c(1,2)]*100))



diff.list <- list(
  RetroPercDiffMin = retro.min.effect,
  RetroPercDiffMax = retro.max.effect,
  RevRetroPercDiffMin = revretro.min.effect,
  RevRetroPercDiffMax = revretro.max.effect,
  JackPercDiffMin = jack.min.effect,
  JackPercDiffMax = jack.max.effect,
  Drop2PercDiff=  cbind(Stock = bm.list$BaseCase[,1],
                        round((bm.list$Drop2[,-c(1,2)] - bm.list$BaseCase[,-c(1,2)])/bm.list$BaseCase[,-c(1,2)]*100))
)

} #end if bm.test


out.list <- list()

if(data.check){ out.list[["DataCheck"]] <- data.check.list}

if(bm.test){out.list[["BM"]] <- bm.list
			out.list[["PercDiff"]] <- diff.list
			}


return(out.list)

}