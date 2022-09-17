
### This file contains all the documentation of datasets ###

#==== smoking - A Smoking Cessation Data ====#
#' @title smoking - A Smoking Cessation Data
#'
#' @aliases smoking
#'
#' @description The original data consist of 223 people enrolled in the study between November 1986
#' and February 1989 from 51 zip codes in the southeastern corner of Minnesota in the
#' United States (Banerjee and Carlin, 2004). In this study, smokers were randomly assigned
#' to one of two treatment groups: smoking intervention (SI) group or usual care (UC) group.
#' The survival time is defined as the time (in years) required for a failed quitter to
#' resume smoking. The people residing in the area with the same zip code form a cluster
#' and may be spatially correlated due to the shared environment.
#'
#' @format Observed covariates include:
#'  \describe{
#'    \item{\code{SexF}}{0 = male, 1 = female.}
#'    \item{\code{Duration}}{duration as smoker in years.}
#'    \item{\code{SI.UC}}{intervention type: 1 = smoking intervention (SI), 0 = usual care (UC).}
#'    \item{\code{F10Cigs}}{the average number of cigarettes smoked per day over the last 10 years (rounded).}
#'    \item{\code{Relapse}}{1 = relapse, 0 = no relapse.}
#'    \item{\code{Timept1}}{the time of study entry.}
#'    \item{\code{Timept2}}{the time of resume smoking.}
#'    \item{\code{Zip}}{51 zip codes in the southeastern corner of Minnesota.}
#'  }
#'
#' @references Banerjee, S. and Carlin,  B. P. (2004) Parametric spatial cure rate models for interval-censored
#' time-to-relapse data. \emph{Biometrics}, \bold{60}: 268-275.
#'
#' @keywords datasets
#'
"smoking"


#==== tonsil - A Smoking Cessation Data ====#
#' @title tonsil - Multi-Center Clinical Trial of Tonsil Carcinoma
#'
#' @aliases tonsil
#'
#' @description A tonsil cancer clinical trial study conducted by the Radiation Therapy Oncology
#' Group in the United States. The survival time is defined as the time (in days) from diagnosis
#' to death. In this study, patients in one institution were randomly assigned to one of two
#' treatment groups: radiation therapy alone or radiation therapy together with a chemotherapeutic
#' agent. A part of the data from the study is available in Kalbfleisch and Prentice (2002).
#'
#' @format A part of the data from the study is available in Kalbfleisch and Prentice (2002),
#' which includes times (in days) from diagnosis to death of 195 patients with squamous cell
#' carcinoma of three sites in the oropharynx between 1968 and 1972 in six participating
#' institutions. Other variables include:
#' \describe{
#'   \item{\code{Inst}}{institution code, from 1 to 6, represents six participating institutions}
#'   \item{\code{Sex}}{1 = male, 2 = female.}
#'   \item{\code{Trt}}{treatment: 1 = standard, 2 = test.}
#'   \item{\code{Grade}}{1 = well differentiated, 2 = moderately differentiated, 3 = poorly differentiated.}
#'   \item{\code{Age}}{in years at time of diagnosis.}
#'   \item{\code{Cond}}{condition: 1 = no disability, 2 = restricted work, 3 = requires assistance with self care, 4 = bed confined.}
#'   \item{\code{Site}}{1 = faucial arch, 2 = tonsillar fossa, 3 = posterior pillar, 4 = pharyngeal tongue, 5 = posterior wall.}
#'   \item{\code{T}}{T staging: 1 = primary tumor measuring 2 cm or less in largest diameter; 2 = primary tumor measuring 2 to 4 cm in largest diameter, minimal infiltration in depth; 3 = primary tumor measuring more than 4 cm; 4 = massive invasive tumor.}
#'   \item{\code{N}}{N staging: 0 = no clinical evidence of node metastases; 1 = single positive node 3 cm or less in diameter, not fixed; 2 = single positive node more than 3 cm in diameter, not fixed; 3 = multiple positive nodes or fixed positive nodes.}
#'   \item{\code{EntryDate}}{Date of entry: Day of year and year.}
#'   \item{\code{Status}}{0 = censored, 1 = dead.}
#'   \item{\code{Time}}{in days from day of diagnosis.}
#' }
#'
#' @references Kalbfleisch, J. D. and Prentice, R. L. (2002) \emph{The Statistical Analysis of Failure Time Data}.
#' John Wiley &  Sons, New York, 2nd edition.
#'
#' @keywords datasets
#'
"tonsil"


#==== bmt - Bone marrow transplantation data ====#
#' @title bmt - Bone marrow transplantation data
#'
#' @aliases bmt
#'
#' @description This multi-center acute leukemia study consists of 137 patients with acute
#' myelocytic leukemia (AML) or acute lymphoblastic leukemia (ALL) aged 7 to 52 from
#' March 1, 1984 to June 30, 1989 at four institutions (Klein and Moeschberger, 2003).
#' The failure time on study is defined at time (in days) to relapse or death.
#'
#' @format The variables represented in the data set are as follows:
#' \describe{
#'   \item{\code{g}}{Disease group: 1 - All, 2 - AML Low Risk, 3 - AML High Risk.}
#'   \item{\code{T1}}{Time to death or on study time.}
#'   \item{\code{T2}}{Disease free survival time (time to relapse, death or end of study).}
#'   \item{\code{d1}}{Death indicator: 1 - Dead, 0 - Alive.}
#'   \item{\code{d2}}{Relapse indicator: 1 - Relapsed, 0 - Disease Free.}
#'   \item{\code{d3}}{Disease free survival indicator: 1 - Dead or Relapsed, 0 - Alive Disease Free.}
#'   \item{\code{TA}}{Time to Acute Graft-Versus-Host Disease.}
#'   \item{\code{da}}{Acute GVHD indicator: 1 - Developed Acute GVHD, 0 - Never Developed Acute GVHD.}
#'   \item{\code{TC}}{Time to Chronic Graft-Versus-Host Disease.}
#'   \item{\code{dc}}{Chronic GVHD Indicator: 1 - Developed Chronic GVHD, 0 - Never Developed Chronic GVHD.}
#'   \item{\code{TP}}{Time to return of platelets to normal levels.}
#'   \item{\code{dp}}{Platelet recovery indicator: 1 - platelets returned to normal, 0 - platelets never returned to normal.}
#'   \item{\code{Z1}}{Patient age in years.}
#'   \item{\code{Z2}}{Donor age in years.}
#'   \item{\code{Z3}}{Patient sex: 1 - Male, 0 - Female.}
#'   \item{\code{Z4}}{Doner sex: 1 - Male, 0 - Female.}
#'   \item{\code{Z5}}{Patient CMV status: 1 - CMV positive, 0 - CMV negative.}
#'   \item{\code{Z6}}{Donor CMV status: 1 - CMV positive, 0 - CMV negative.}
#'   \item{\code{Z7}}{Waiting time to transplant in days.}
#'   \item{\code{Z8}}{FAB: 1 - FAB Grade 4 or 5 and AML, 0 - otherwise.}
#'   \item{\code{Z9}}{Hospital: 1 - The Ohio State University, 2 - Alferd , 3 - St. Vincent, 4 - Hahnemann.}
#'   \item{\code{Z10}}{MTX: used as a Graft-Versus-Host-Prophylactic 1 - Yes, 0 - No.}
#' }
#'
#' @references Klein, J. P. and Moeschberger, M. L. (2003) \emph{Survival Analysis: Techniques for Censored and Truncated Data}.
#' Springer, New York, 2nd edition.
#'
#' @keywords datasets
#'
"bmt"



