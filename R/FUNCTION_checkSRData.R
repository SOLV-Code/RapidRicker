#' checkSRData
#'
#' This function calculates a set of diagnostics for the spawner-recruit data. Some apply to the whole series 
#'  (e.g. the contrast is spawner estimates), while others flag individual observations
#'  (e.g. R/S above user-specified plausible upper bound, pointing to a potential data error in either R or S).
#' @param sr_obj  a data frame with  Year, Spn, Rec, and optionally SpnExp and RecAge#. (Data for 1 Stock!). Other variables can be there but are not used.
#' @param flags a data frame with Label, Lower, Upper. Lower and Upper define the triggers for values to be flagged. Labels have to match the criteria used in the function(LIST). The built-in object flags_default has is a template. 
#' @param trace if TRUE, print various intermediate details to the screen
#' @keywords contrast, outliers
#' @export
#' @examples
#' data.chk <- checkSRData(SR_Sample[SR_Sample$Stock == "Stock1",])
#' print(data.chk)


checkSRData <- function(sr_obj,flags = NULL,trace = FALSE){



req.inputs <- c("Year", "Spn", "Rec") 

if(!all(req.inputs %in% names(sr_obj))){warning("in call to checkSRData(): One of the required columns is missing: Year, Spn, Rec"); stop()}


if(is.null(flags)){ flags <- flags_default }  # use built in data object unless user specifies different


# Low contrast in Spn Data -------------------------------------
# calc: max(spn)/min(spn)

if(trace){print("Starting Contr")}

flag.label <- "Contr"
metric.in <- max(sr_obj$Spn,na.rm=TRUE)/min(sr_obj$Spn,na.rm=TRUE)

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                       metric.val = metric.in,
                       round.val = 2)

flags[flags$Label == flag.label, "MetricVal"] <- flag.tmp$MetricVal
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged




# Insufficient SR data -------------------------------------
# calc: num(BY) with spn and rec

if(trace){print("Starting numObs")}

sr.idx <- !is.na(sr_obj$Spn) & !is.na(sr_obj$Rec)
metric.in  <- sum(sr.idx)
flag.label <- "NumObs"

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                      metric.val = metric.in,
                      round.val = 2)

flags[flags$Label == flag.label, "MetricVal"] <- flag.tmp$MetricVal
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged




# Missing large Spn -------------------------------------
# calc: Max(spn)/(Max(spn in SR)

if(trace){print("Starting LgSpn")}

metric.in<- max(sr_obj$Spn,na.rm = TRUE) / max(sr_obj$Spn[sr.idx],na.rm = TRUE)
flag.label <- "LgSpn"

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                      metric.val = metric.in,
                      round.val = 2)

flags[flags$Label == flag.label, "MetricVal"] <- flag.tmp$MetricVal
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged



# Missing large Rec -------------------------------------
# calc: Max(Rec)/(Max(Rec in SR)

if(trace){print("Starting LgRec")}

metric.in <- max(sr_obj$Rec,na.rm = TRUE) / max(sr_obj$Rec[sr.idx],na.rm = TRUE)
flag.label <- "LgRec"

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                      metric.val = metric.in,
                      round.val = 2)

flags[flags$Label == flag.label, "MetricVal"] <- flag.tmp$MetricVal
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged




# Missing small Spn -------------------------------------
# calc: Min(spn)/(Min(spn in SR)

if(trace){print("Starting SmSpn")}

metric.in<- min(sr_obj$Spn,na.rm = TRUE) / min(sr_obj$Spn[sr.idx],na.rm = TRUE)
flag.label <- "SmSpn"

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                      metric.val = metric.in,
                      round.val = 2)

flags[flags$Label == flag.label, "MetricVal"] <- flag.tmp$MetricVal
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged


# Missing small Rec -------------------------------------
# calc: Min(rec)/(Min(rec in SR)

if(trace){print("Starting SmRec")}

metric.in<- min(sr_obj$Rec,na.rm = TRUE) / min(sr_obj$Rec[sr.idx],na.rm = TRUE)
flag.label <- "SmRec"

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                      metric.val = metric.in,
                      round.val = 2)

flags[flags$Label == flag.label, "MetricVal"] <- flag.tmp$MetricVal
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged




# Large Spn Expansion Factor (median) --------------------------------------------
# calc: Median(ExpSpn)

if(trace){print("Starting LgExp")}

flag.label <- "LgExp"

if("SpnExp" %in% names(sr_obj)){

metric.in <- median(sr_obj$SpnExp, na.rm = TRUE)


flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                      metric.val = metric.in,
                      round.val = 2)

flags[flags$Label == flag.label, "MetricVal"] <- flag.tmp$MetricVal
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged
}

if(!("SpnExp" %in% names(sr_obj))){
  flags[flags$Label == flag.label, "MetricVal"] <- NA
  flags[flags$Label == flag.label, "Flagged"] <- NA
}




# Unstable age comp --------------------------------------------
# calc: MaxDiffByAge

if(FALSE){ # excluded, until resolve this issue: https://github.com/SOLV-Code/RapidRicker/issues/61

if(trace){print("Starting VarAge")}

flag.label <- "VarAge"

if(sum(grepl("RecAge",names(sr_obj)))>=2) { # don only of have at least 2 age classes in the data


  age.sub <- sr_obj[,grepl("RecAge", names(sr_obj))]

  metric.pre <- lapply(age.sub,range.calc)

  metric.in <- max(unlist(metric.pre))

  flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                        flag.upper = flags[flags$Label == flag.label,"Upper"],
                        metric.val = metric.in,
                        round.val = 2)

  flags[flags$Label == flag.label, "MetricVal"] <- flag.tmp$MetricVal
  flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged
}

if(sum(grepl("RecAge",names(sr_obj)))<2){
  flags[flags$Label == flag.label, "MetricVal"] <- NA
  flags[flags$Label == flag.label, "Flagged"] <- NA
}

}


# Unusual Spn obs --------------------------------------------
# calc: Spn / Med(Spn)


if(trace){print("Starting OddSpn")}

flag.label <- "OddSpn"
metric.obs <- sr_obj$Spn / median(sr_obj$Spn, na.rm = TRUE)

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                        flag.upper = flags[flags$Label == flag.label,"Upper"],
                        metric.val = metric.obs,
                        round.val = 2)
flags[flags$Label == flag.label, "MetricVal"] <- NA
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged
flags[flags$Label == flag.label, "NumFlagged"] <-sum(flag.tmp$FlaggedObs,na.rm=TRUE)

sr_obj[[paste0(flag.label,"Val")]] <- flag.tmp$MetricVal
sr_obj[[paste0(flag.label,"Flag")]] <- flag.tmp$FlaggedObs



# Unusual Rec obs --------------------------------------------
# calc: Rec / Med(Rec)

if(trace){print("Starting OddRec")}

flag.label <- "OddRec"
metric.obs <- sr_obj$Rec/ median(sr_obj$Rec, na.rm = TRUE)

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                      metric.val = metric.obs,
                      round.val = 2)
flags[flags$Label == flag.label, "MetricVal"] <- NA
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged
flags[flags$Label == flag.label, "NumFlagged"] <-sum(flag.tmp$FlaggedObs,na.rm=TRUE)

sr_obj[[paste0(flag.label,"Val")]] <- flag.tmp$MetricVal
sr_obj[[paste0(flag.label,"Flag")]] <- flag.tmp$FlaggedObs



# Unusual Prod obs --------------------------------------------
# calc: R/S (or scale by median R/S? NEED TO DISCUSS

if(trace){print("Starting OddProd")}

flag.label <- "OddProd"
metric.obs <- (sr_obj$Rec/sr_obj$Spn) #/ median(sr_obj$Rec/sr_obj$Spn, na.rm = TRUE)

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                      metric.val = metric.obs,
                      round.val = 2)
flags[flags$Label == flag.label, "MetricVal"] <- NA
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged
flags[flags$Label == flag.label, "NumFlagged"] <-sum(flag.tmp$FlaggedObs,na.rm=TRUE)

sr_obj[[paste0(flag.label,"Val")]] <- flag.tmp$MetricVal
sr_obj[[paste0(flag.label,"Flag")]] <- flag.tmp$FlaggedObs


# Unusual Spn Expansions --------------------------------------------
# calc: SpnExp / Med(SpnExp)

if(trace){print("Starting OddExp")}

flag.label <- "OddExp"

if("SpnExp" %in% names(sr_obj)){

metric.obs <- sr_obj$SpnExp   / median(sr_obj$SpnExp, na.rm = TRUE)

flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                      flag.upper = flags[flags$Label == flag.label,"Upper"],
                      metric.val = metric.obs,
                      round.val = 2)
flags[flags$Label == flag.label, "MetricVal"] <- NA
flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged
flags[flags$Label == flag.label, "NumFlagged"] <-sum(flag.tmp$FlaggedObs,na.rm=TRUE)

sr_obj[[paste0(flag.label,"Val")]] <- flag.tmp$MetricVal
sr_obj[[paste0(flag.label,"Flag")]] <- flag.tmp$FlaggedObs
} # end if have SpnExp column



if(!("SpnExp" %in% names(sr_obj))){

flags[flags$Label == flag.label, "MetricVal"] <- NA
flags[flags$Label == flag.label, "Flagged"] <- NA
flags[flags$Label == flag.label, "NumFlagged"] <- NA

sr_obj[[paste0(flag.label,"Val")]] <- NA
sr_obj[[paste0(flag.label,"Flag")]] <- NA


}





# Unusual Age Comp --------------------------------------------
# calc: recAge / Med(RecAge) for main age comp
#


if(FALSE){ # excluded, until resolve this issue: https://github.com/SOLV-Code/RapidRicker/issues/61

if(trace){print("Starting OddAge")}

flag.label <- "OddAge"

if(sum(grepl("RecAge",names(sr_obj)))>=2) { # do only of have at least 2 age classes in the data

  age.sub <- sr_obj[,grepl("RecAge", names(sr_obj))]
  metric.pre <- lapply(age.sub,function(x){median(x,na.rm =TRUE)})
  age.main <- names(unlist(metric.pre))[unlist(metric.pre) == max(unlist(metric.pre))]

  metric.in <- sr_obj[[age.main]] / metric.pre[[age.main]]

  flag.tmp <- flag.calc(flag.lower = flags[flags$Label == flag.label,"Lower"],
                        flag.upper = flags[flags$Label == flag.label,"Upper"],
                        metric.val = metric.in,
                        round.val = 2)

  flags[flags$Label == flag.label, "MetricVal"] <- NA
  flags[flags$Label == flag.label, "Flagged"] <- flag.tmp$Flagged
  flags[flags$Label == flag.label, "NumFlagged"] <-sum(flag.tmp$FlaggedObs,na.rm=TRUE)

  sr_obj[[paste0(flag.label,"Val")]] <- flag.tmp$MetricVal
  sr_obj[[paste0(flag.label,"Flag")]] <- flag.tmp$FlaggedObs
}

if(sum(grepl("RecAge",names(sr_obj)))<2){
  flags[flags$Label == flag.label, "MetricVal"] <- NA
  flags[flags$Label == flag.label, "Flagged"] <- NA
  flags[flags$Label == flag.label, "NumFlagged"] <- NA
  sr_obj[[paste0(flag.label,"Val")]] <- NA
  sr_obj[[paste0(flag.label,"Flag")]] <- NA
}


}


#-----------------------------------------

main.cols <- flags$Scope == "Obs" & !grepl("Age",flags$Label)

sr_obj[["NumFlagsMain"]]  <- rowSums(sr_obj[,paste0(flags$Label[main.cols],"Flag")])
sr_obj[["NumFlagsAll"]]  <- rowSums(sr_obj[,paste0(flags$Label[flags$Scope == "Obs"],"Flag")])


list.out <- list(Summary = flags, Data = sr_obj)

} # end checkSRData() function




flag.calc <- function(flag.lower, flag.upper,metric.val,round.val = NULL){
# internal subroutine
# metric.val can be a single value or a vector

  # based on https://stackoverflow.com/questions/16822426/dealing-with-true-false-na-and-nan
  flag.out <- any(c(metric.val < flag.lower, metric.val > flag.upper)  %in% TRUE)
  #print(metric.val < flag.lower)
  #print(metric.val > flag.upper)


  if(length(metric.val)>1){flag.vec <- ((metric.val < flag.lower) %in% TRUE) |   ((metric.val > flag.upper) %in% TRUE)  }
  if(length(metric.val)==1){flag.vec <- NA}


  if(!is.null(round.val)){ val.out <- round(metric.val,2) }
  if(is.null(round.val)){ val.out <- metric.val }

  out.list <- list(MetricVal = val.out, Flagged = flag.out, FlaggedObs = flag.vec)
  return(out.list)
}


range.calc <- function(x){
# internal subroutine
  range.out <- max(x,na.rm = TRUE) - min(x,na.rm = TRUE)
  return(range.out)
}


div.by.median <- function(x){
# internal subroutine
  vec.out <- x / median(x,na.rm=TRUE)
}