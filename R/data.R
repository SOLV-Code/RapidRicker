#' Sample Spawner Recruit Data
#'
#' A dataset containing spawners, recruits, spawner estimate expansions, and age comp by brood year for mutiple stocks.
#'
#' @format A data frame with  5 required and some optional variables:
#' \describe{
#'   \item{Stock}{ stock name}
#'   \item{Year}{brood year}
#'   \item{Spn}{spawner estimate}
#'   \item{Rec}{recruitment estimate}
#'   \item{logRpS}{calculated value of log(Rec/Spn), usually with mutate(logRpS = log(Rec/Spn))}
#'   \item{SpnExp}{OPTIONAL: expansion factor applied to get the spawner estimate}
#'   \item{RecAgeNumber}{proportion each age class contributes to the total recruits (e.g. RecAge3, RecAge4, RecAge5}
#'   \item{}{}
#'   ...
#' }
#' @source Dummy Data
"SR_Sample"


#' Upper and Lower thresholds for flagging SR data issues
#'
#' A dataset containing labels, thresholds and decriptions for various data qulity metrics. INCLUDE LIST!
#'
#' @format A data frame with  3 required variables and optional context information (e.g. Reference):
#' \describe{
#'   \item{Label}{Labels identifyin the criteria. Have to match the criteria used in the function.}
#'   \item{Lower}{Any value lower than this is flagged.}
#'   \item{Upper}{Any value higher than this is flagged.}
#'   \item{}{}
#'   ...
#' }
#' @source Include REF
"flags_default"


#' Spawner Recruit Test Data
#'
#' A dataset containing spawners and  recruits by brood year for 2 stocks.
#'
#' @format A data frame with  5 required variables:
#' \describe{
#'   \item{Stock}{ stock name}
#'   \item{Year}{brood year}
#'   \item{Spn}{spawner estimate}
#'   \item{Rec}{recruitment estimate}
#'   \item{logRpS}{calculated value of log(Rec/Spn), usually with mutate(logRpS = log(Rec/Spn))}
#'   \item{}{}
#'   ...
#' }
#' @source Dummy Data
"gls_test_data"



