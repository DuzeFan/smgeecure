

#=================================#
# The overall package description #
#=================================#

#' smgeecure: Mixture Cure Models with Auxiliary Subgroup Survival Probabilities.
#'
#' This package features the marginal semi-parametric proportional hazards/accelerated failure time mixture cure models for analyzing clustered survival data with a possible cure fraction.
#'    The underlying methods are based on the paper titled "Marginal Accelerated Failure Time Cure Model for clustered survival data with long-term survivors",
#'    which has been submitted to the journal: Statistical Methods in Medical Research.
#'
#' @importFrom methods is
#' @importFrom stats cov model.frame model.matrix model.response optim pchisq pnorm rbinom rexp runif sd uniroot
#' @docType package
#' @name smgeecure
NULL



#========================================================#
# Fit the semiparametric marginal mixture cure model ####
#========================================================#

#==== The main function for smc aft/ph mixture cure fitting model ====#
#' @title Marginal semi-parametric mixture cure model with generalized estimating equations
#'
#' @description Fit the semiparametric marginal proportional hazards mixture cure (PHMC) model or
#' the semiparametric accelerated failure time mixture cure (AFTMC) model with the generalized
#' estimating equations (GEE). GEE approach is generalized through the expectation-solution (ES) algorithm
#' to account for the correlation among the cure statuses and the dependence among the failure
#' times of uncured patients to improve the estimation efficiency.
#'
#' @aliases smgeecure
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}.
#'   The \code{response} is a \code{Surv} object with right censoring. It is used to specify the covariate effects on the failure time of uncured subjects. See the documentation for \code{survreg} and \code{Surv} in package \code{survival} for details.
#'   The expression to the right of the "~" specifies the effect of covariates on the failure time of uncured patients.
#' @param cureform indicator function a formula expression, of the form \code{cureform ~ predictors}. It is used to specify the effects of covariates on the cure rate. A covariate may be used in both \code{formula} and \code{cureform}.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param data a data frame in which to interpret the variables named in the \code{formula} and the \code{cureform}.
#' @param model specifies the model using in latency part. It can be \code{ph} which represents
#' the proportional hazards model, or \code{aft} which represents the accelerated failure time model.
#' @param corstr a character string specifying the correlation structure. Currently, the following three working correlation structures are available:
#' \code{independence} for the independent working correlation structure,
#' \code{exchangeable} for the exchangeable working correlation structure and
#' \code{ar1} for the first-order auto-regressive structure.
#' @param Var a logical value. If it is TRUE, the program returns standard error estimates for beta and gamma.
#' Otherwise, no standard errors are computed and only coefficient estimates are returned.
#' By default, \code{Var = TRUE}. Note that for both \code{model="ph"} and \code{"aft"}, the returned
#' standard error estimators are based on the bootstrap method.
#' @param nboot specifies the number of bootstrap samplings. The default value is 100.
#' @param stdz a logical value. If it is \code{TRUE}, all the covariates are standardized. Otherwise, non-transformed covariates are used.
#' @param esmax specifies the maximum iteration number. If the convergence criterion is not met,
#' the ES iteration will be stopped after \code{esmax} iterations and the estimates will be based on the
#' last ES iteration. The default \code{esmax = 100}.
#' @param eps specifies the relative tolerance, and is set to \code{1e-06} as the default.
#' The iterations are considered to be converged when the maximum relative change in the
#' parameters and likelihood estimates between iterations is less than the value specified.
#'
#' @details Semiparametric marginal PHMC and AFTMC models are considered in this function.
#' For cure rate, a logistic regression model is employed and the probability of being cured is given
#' by \code{(1+exp(gamZ))^{(-1)}}.
#' A covariate can be used either in formula or in cureform or in both. The model parameters are
#' estimated by the expectation-solution (ES) algorithm and the standard error estimates are obtained
#' using bootstrap.
#'
#' @return An object of class \code{smgeecure} is returned. It can be examined by \code{print.smgeecure()}.
#'
#' @references Cai, C., Zou, Y., Peng, Y., and Zhang, J. (2012). smcure: An r-package for estimating semiparametric
#' mixture cure models. Computer methods and programs in biomedicine, 108(3):1255.
#' @references Zhang, J. and Peng, Y. (2007). A new estimation method for the semiparametric accelerated
#' failure time mixture cure model. Statistics in medicine, 26(16):3157–3171.
#' @references Niu, Y., Song, L., Liu, Y, and Peng, Y. (2018) Modeling clustered long-term survivors using marginal
#' mixture cure model. Biometrical Journal, doi: 10.1002/bjmj.201700114.
#'
#' @examples
#' \donttest{
#' ### Be patient, the following examples may take several minutes on a faster computer.
#'
#' ### library
#' library(smgeecure)
#' library(geecure)
#'
#' #=======================================================================#
#' ### Example 1. Fit the marginal semiparametric PHMC model for
#' ###             the tonsil data.
#'
#' # ---- prepare the data
#' data(tonsil)
#' tonsil <- tonsil[-c(141,136,159),]
#' tonsil$Sex <- ifelse(tonsil$Sex == 1, 0, 1) # 1="Female"
#' tonsil$Cond <- ifelse(tonsil$Cond == 1, 0, 1) # 0=no disability
#' tonsil$T <- ifelse(tonsil$T < 4, 0, 1)
#' tonsil$Grade2 <- ifelse(tonsil$Grade==2,1,0)
#' tonsil$Grade3 <- ifelse(tonsil$Grade==3,1,0)
#' table(tonsil$Inst)
#'
#' # ---- plot the KM survival curve (overall)
#' pdf("Figure_Tonsil_KM.pdf",width=7,height=6)
#' plot(survival::survfit(survival::Surv(Time, Status) ~ 1, data = tonsil),
#'      conf.int = T, mark.time = TRUE, lwd=1.5,
#'      ylab = "Survival Probability", xlab = "Survival Time (in Days)",
#'      xlim = c(0,2000), ylim = c(0,1),
#'      cex.lab = 1,
#'      xaxt = "n", yaxt = "n"
#' )
#' axis(1,seq(0,2000,500),seq(0,2000,500),cex.axis = 1)
#' axis(2,seq(0,1,0.2),seq(0,1,0.2),cex.axis = 1)
#' dev.off()
#'
#' # ---- plot the KM survival curve (stratified)
#' pdf("Figure_Tonsil_KM_SexTumorsize.pdf",width=7,height=6)
#' plot(survival::survfit(survival::Surv(Time, Status) ~ Sex + T, data = tonsil),
#'      conf.int = F,  mark.time = TRUE, lwd=1.5,
#'      ylab = "Survival Probability", xlab = "Survival Time (in Days)",
#'      xlim = c(0,2000), ylim = c(0,1),
#'      lty = c(1,2,3,4), col = c(1,2,3,4),
#'      xaxt = "n", yaxt = "n", cex.lab=1)
#' axis(1,seq(0,2000,400),seq(0,2000,400),cex.axis=1)
#' axis(2,seq(0,1,0.2),seq(0,1,0.2),cex.axis=1)
#' legend("topright", c("Massive/Female","Massive/Male","Non-massive/Female","Non-massive/Male"),
#'        cex = 1, col = c(1,2,3,4), lty = c(1,2,3,4))
#' dev.off()
#'
#' # ---- plot the log-cumulative survival curve (check ph assumption)
#' #         Plot BY Sex
#' # estimates for 1
#' KM.fit.1 <- survival::survfit(survival::Surv(Time, Status) ~ 1,
#'                               data = tonsil, subset = (Sex==1))
#' Infty.1 <- min(KM.fit.1$surv)
#' surv.1.Latency <- (KM.fit.1$surv-Infty.1) / (1-Infty.1)
#' time.1.Latency <- KM.fit.1$time
#' cumu.1.Latency.Nozero <- -log(surv.1.Latency[surv.1.Latency!=0])
#' time.1.Latency.Nozero <- KM.fit.1$time[surv.1.Latency!=0]
#' # estimates for 0
#' KM.fit.0 <- survival::survfit(survival::Surv(Time, Status) ~ 1,
#'                               data = tonsil, subset = (Sex==0))
#' Infty.0 <- min(KM.fit.0$surv)
#' surv.0.Latency <- (KM.fit.0$surv-Infty.0) / (1-Infty.0)
#' time.0.Latency <- KM.fit.0$time
#' cumu.0.Latency.Nozero <- -log(surv.0.Latency[surv.0.Latency!=0])
#' time.0.Latency.Nozero <- KM.fit.0$time[surv.0.Latency!=0]
#' # plot of log-cumulative function
#' pdf("Figure_Tonsil_LogCum_Sex.pdf",width=7,height=6)
#' plot(log(cumu.1.Latency.Nozero)~log(time.1.Latency.Nozero),type="s",lty=1,col="red",
#'      xlim = c(0,8), ylim=c(-5,1.5), lwd=1.5,
#'      ylab = "Logarithm of the Cumulative Hazard Function",
#'      xlab = "Logarithm of the Survival Time",
#'      xaxt = "n",
#'      cex.lab=1
#' )
#' lines(log(cumu.0.Latency.Nozero)~log(time.0.Latency.Nozero),type="s",lty=2,col="black",lwd=1.5)
#' axis(1,seq(0,8,2),seq(0,8,2),cex.axis = 1)
#' legend("bottomright", c("Female","Male"), cex = 1,
#'        col = c("red","black"), lty = c(1,2))
#' dev.off()
#'
#' # ---- plot the log-cumulative survival curve (check ph assumption)
#' #         Plot BY Sex
#' #         (Just use uncensored patients)
#' UncuredPatients <- tonsil$Status==1
#' pdf("Figure_Tonsil_LogCum_Uncen.pdf",width=7,height=6)
#' xlim_lower <- 0; xlim_upper <- 1800; xlim_step  <- 300
#' ylim_lower <- -5; ylim_upper <- 2; ylim_step  <- 1
#' mylty <- c(1,2); mycol <- c("black","red")
#' plot(survival::survfit(survival::Surv(Time, Status) ~ Sex,
#'                        data = tonsil, subset=UncuredPatients,
#'                        stype=1),
#'      conf.int = F,  mark.time = FALSE, cumhaz = TRUE, log="xy",
#'      ylab = "Logarithm of the Cumulative Hazard Function",
#'      xlab = "Logarithm of the Survival Time (in Days)",
#'      # xlim = c(xlim_lower,xlim_upper), ylim = c(ylim_lower,ylim_upper),
#'      lty = mylty, col = mycol, lwd=1.5,
#'      xaxt = "n", yaxt = "n")
#' axis(1,c(100,300,600,1500),c(100,300,600,1500))
#' axis(2,exp(seq(ylim_lower,ylim_upper,ylim_step)),seq(ylim_lower,ylim_upper,ylim_step))
#' legend("bottomright", c("male","female"), cex = 1.3, #bty="n",
#'        col = mycol, lty = mylty)
#' dev.off()
#'
#' # ---- Fit the data using marginal semi-parametric marginal AFTMC model
#'
#' # fit marginal semi-parametric marginal AFTMC model (exchangeable correlation)
#' set.seed(521)
#' tonsil.aft.gee.ex <- smgeecure(
#'   formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T,
#'   cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst,
#'   data = tonsil, model = "aft", corstr = "exchangeable", Var = T, nboot = 100
#' )
#' print.smgeecure(tonsil.aft.gee.ex)
#'
#' # fit marginal semi-parametric AFTMC model (independence correlation)
#' set.seed(521)
#' tonsil.aft.gee.ind <- smgeecure(
#'   formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T,
#'   cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst,
#'   data = tonsil, model = "aft", corstr = "independence", Var = T, nboot = 100
#' )
#' print.smgeecure(tonsil.aft.gee.ind)
#'
#' # output the results
#' write.csv(cbind(tonsil.aft.gee.ex$incidence,rbind(NA,tonsil.aft.gee.ex$latency)),
#'           "tonsil.aft.gee.ex.csv")
#' write.csv(cbind(tonsil.aft.gee.ind$incidence,rbind(NA,tonsil.aft.gee.ind$latency)),
#'           "tonsil.aft.gee.ind.csv")
#'
#' # ---- Fit the data using marginal semi-parametric marginal PHMC model
#'
#' # fit marginal semi-parametric PHMC model (exchangeable correlation)
#' set.seed(1)
#' tonsil.ph.gee.ex <- smgeecure(
#'   formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T,
#'   cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst,
#'   data = tonsil, model = "ph", corstr = "exchangeable", Var = T,nboot = 100
#' )
#' print.smgeecure(tonsil.ph.gee.ex)
#'
#' # fit marginal semi-parametric PHMC model (independence correlation)
#' set.seed(1)
#' tonsil.ph.gee.ind <- smgeecure(
#'   formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T,
#'   cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst,
#'   data = tonsil, model = "ph", corstr = "independence", Var = F, nboot = 100
#' )
#' print.smgeecure(tonsil.ph.gee.ind)
#'
#' # output the results
#' write.csv(cbind(tonsil.ph.gee.ex$incidence,rbind(NA,tonsil.ph.gee.ex$latency)),
#'           "tonsil.ph.gee.ex.csv")
#' write.csv(cbind(tonsil.ph.gee.ind$incidence,rbind(NA,tonsil.ph.gee.ind$latency)),
#'           "tonsil.ph.gee.ind.csv")
#'
#'
#' #=======================================================================#
#' ### Example 2. Fit the semiparametric marginal PHMC model for
#' ###             the bone marrow transplantation data.
#'
#' # ---- prepare the data
#' data(bmt)
#' bmt$g <- factor(bmt$g, label = c("ALL", "AML low risk","AML high risk"))
#' bmt$Z8 <- factor(bmt$Z8, label = c("Otherwise", "FAB"))
#'
#' # ---- plot the KM survival curve (overall and stratified)
#' pdf("Figure_bmt_KM.pdf",width=7,height=6)
#' plot(survival::survfit(survival::Surv(T2, d3) ~ 1, data = bmt),
#'      conf.int = T, mark.time = TRUE, lwd=1.5,
#'      ylab = "Survival Probability", xlab = "Survival Time (in Days)",
#'      xlim = c(0,2800), ylim = c(0,1),
#'      xaxt = "n", yaxt = "n",
#'      cex.lab = 1
#' )
#' axis(1,seq(0,2800,400),seq(0,2800,400),cex.axis = 1)
#' axis(2,seq(0,1,0.2),seq(0,1,0.2),cex.axis = 1)
#' dev.off()
#'
#' # ---- plot the KM survival curve (overall and stratified)
#' pdf("Figure_bmt_KM_g.pdf",width=7,height=6)
#' plot(survival::survfit(survival::Surv(T2, d3) ~ factor(g), data = bmt),
#'      conf.int = F,  mark.time = TRUE, lwd=1.5,
#'      ylab = "Survival Probability", xlab = "Survival Time (in Days)",
#'      xlim = c(0,2800), ylim = c(0,1),
#'      lty = c(1,2,3), col = c(1,2,3),
#'      xaxt = "n", yaxt = "n", cex.lab=1)
#' axis(1,seq(0,2800,400),seq(0,2800,400),cex.axis=1)
#' axis(2,seq(0,1,0.2),seq(0,1,0.2),cex.axis=1)
#' legend("topright", c("ALL","AML low risk","AML high risk"),
#'        cex = 1, col = c(1,2,3), lty = c(1,2,3))
#' dev.off()
#'
#'
#' # ---- plot the log-cumulative survival curve (check ph assumption)
#' #         Plot by g's two dummy factors
#' pdf("Figure_bmt_LogCum_g.pdf",width=10,height=5)
#' par(mfrow = c(1,2))
#' groups.names <- c("AML low risk","AML high risk")
#' groups <- cbind(bmt$g==groups.names[1],bmt$g==groups.names[2])
#' for(igroup in 1:ncol(groups)){
#'   # estimates for 1
#'   KM.fit.1 <- survival::survfit(survival::Surv(T2, d3) ~ 1,
#'                                 data = bmt, subset = groups[,igroup])
#'   Infty.1 <- min(KM.fit.1$surv)
#'   surv.1.Latency <- (KM.fit.1$surv-Infty.1) / (1-Infty.1)
#'   time.1.Latency <- KM.fit.1$time
#'   cumu.1.Latency.Nozero <- -log(surv.1.Latency[surv.1.Latency!=0])
#'   time.1.Latency.Nozero <- KM.fit.1$time[surv.1.Latency!=0]
#'   # estimates for 0
#'   KM.fit.0 <- survival::survfit(survival::Surv(T2, d3) ~ 1,
#'                                 data = bmt, subset = !groups[,igroup])
#'   Infty.0 <- min(KM.fit.0$surv)
#'   surv.0.Latency <- (KM.fit.0$surv-Infty.0) / (1-Infty.0)
#'   time.0.Latency <- KM.fit.0$time
#'   cumu.0.Latency.Nozero <- -log(surv.0.Latency[surv.0.Latency!=0])
#'   time.0.Latency.Nozero <- KM.fit.0$time[surv.0.Latency!=0]
#'   # plot of log-cumulative function
#'   plot(log(cumu.1.Latency.Nozero)~log(time.1.Latency.Nozero),type="s",lty=1,col="red",
#'        xlim = c(0,8), ylim=c(-4.5,1.5), lwd=1.5,
#'        ylab = "Logarithm of the Cumulative Hazard Function",
#'        xlab = "Logarithm of the Survival Time",
#'        main=groups.names[igroup],
#'        xaxt = "n",
#'        cex.lab = 1,cex.main=1
#'   )
#'   lines(log(cumu.0.Latency.Nozero)~log(time.0.Latency.Nozero),type="s",lty=2,col="black",lwd=1.5)
#'   axis(1,seq(0,8,2),seq(0,8,2),cex.axis = 1)
#'   legend("bottomright", c(groups.names[igroup],"Other"), cex = 1,
#'          col = c("red","black"), lty = c(1,2))
#' }
#' dev.off()
#'
#'
#' # ---- Fit the data using marginal semi-parametric marginal PHMC model
#'
#' # fit marginal semi-parametric PHMC model (exchangeable correlation)
#' set.seed(1)
#' bmt.ph.gee.ex <- smgeecure(
#'   formula = Surv(T2, d3) ~ factor(g) + Z8, cureform = ~ factor(g) + Z8,
#'   id = bmt$Z9, data = bmt, model = "ph", corstr = "exchangeable",
#'   Var = T,nboot = 100, esmax = 100, eps = 1e-06
#' )
#' print.smgeecure(bmt.ph.gee.ex)
#'
#' # fit marginal semi-parametric PHMC model (independence correlation)
#' set.seed(1)
#' bmt.ph.gee.ind <- smgeecure(
#'   formula = Surv(T2, d3) ~ factor(g) + Z8, cureform = ~ factor(g) + Z8,
#'   id = bmt$Z9, data = bmt, model = "ph", corstr = "independence",
#'   Var = T, nboot = 100, esmax = 100, eps = 1e-06
#' )
#' print.smgeecure(bmt.ph.gee.ind)
#'
#' # output the results
#' write.csv(cbind(bmt.ph.gee.ex$incidence,rbind(NA,bmt.ph.gee.ex$latency)),
#'           "bmt.ph.gee.ex.csv")
#' write.csv(cbind(bmt.ph.gee.ind$incidence,rbind(NA,bmt.ph.gee.ind$latency)),
#'           "bmt.ph.gee.ind.csv")
#'
#' # ---- Fit the data using marginal semi-parametric marginal AFTMC model
#'
#' # fit marginal semi-parametric marginal AFTMC model (exchangeable correlation)
#' set.seed(1)
#' bmt.aft.gee.ex <- smgeecure(
#'   formula = Surv(T2, d3) ~ factor(g) + Z8, cureform = ~ factor(g) + Z8,
#'   id = bmt$Z9, data = bmt, model = "aft", corstr = "exchangeable",
#'   Var = T,nboot = 100
#' )
#' print.smgeecure(bmt.aft.gee.ex)
#'
#' # fit marginal semi-parametric AFTMC model (independence correlation)
#' set.seed(1)
#' bmt.aft.gee.ind <- smgeecure(
#'   formula = Surv(T2, d3) ~ factor(g) + Z8, cureform = ~ factor(g) + Z8,
#'   id = bmt$Z9, data = bmt, model = "aft", corstr = "independence",
#'   Var = T,nboot = 100
#' )
#' print.smgeecure(bmt.aft.gee.ind)
#'
#' # output the results
#' write.csv(cbind(bmt.aft.gee.ex$incidence,rbind(NA,bmt.aft.gee.ex$latency)),
#'           "bmt.aft.gee.ex.csv")
#' write.csv(cbind(bmt.aft.gee.ind$incidence,rbind(NA,bmt.aft.gee.ind$latency)),
#'           "bmt.aft.gee.ind.csv")
#' }
#'
#' @export smgeecure
smgeecure <- function(
    formula, cureform, id, data, model = c("aft","ph"),
    corstr = c("independence", "exchangeable","ar1"),
    Var = TRUE, nboot = 100, stdz = TRUE, esmax = 20, eps = 1e-04
){

  ## basic
  call <- match.call()
  Y <- model.response(model.frame(formula,data))
  if (!inherits(Y, "Surv")){stop("Response must be a survival object")}
  uid <- sort(unique(id))

  ## fit the model in various situations
  if(model=="ph"){

    res <- geecure(formula=formula,cureform=cureform,data=data,id=id,
                   model="semi",corstr=corstr,Var=FALSE,boots=Var,nboot=nboot,
                   stdz=stdz,esmax=esmax,eps=eps)
    class(res) <- "list"


  }else if(model=="aft"){

    ## prepare data - rearrange our subjects according to ids
    uid <- sort(unique(id))
    idx.rearrange <- unlist(lapply(uid,function(iid){which(id == iid)}))
    id.rearrange <- unlist(lapply(1:length(uid),function(ii){rep(ii,sum(id == uid[ii]))}))
    data <- data[idx.rearrange,,drop=F]
    id <- id.rearrange
    rownames(data) <- NULL

    # prepare - necesary vectors and matrices + extract vars from formulas
    mf.latency <- model.frame(formula,data)
    mf.incidence <- model.frame(cureform,data)
    X <-  model.matrix(attr(mf.latency, "terms"), mf.latency)[,-1,drop=F]
    Z <-  model.matrix(attr(mf.incidence, "terms"), mf.incidence)[,-1,drop=F]
    vars.name.response <- all.vars(formula)[c(1,2)]
    vars.name.latency <- sapply(colnames(X),function(x){for(xx in c(" ","(",")",":")){x <- gsub(xx,"_",x,fixed=T)};x},USE.NAMES=F)
    vars.name.incidence <- sapply(colnames(Z),function(x){for(xx in c(" ","(",")",":")){x <- gsub(xx,"_",x,fixed=T)};x},USE.NAMES=F)
    vars.name.covariate <- union(vars.name.latency,vars.name.incidence)
    colnames(X) <- vars.name.latency; colnames(Z) <- vars.name.incidence
    tobs  <- data[,vars.name.response[1]]
    delta <- data[,vars.name.response[2]]

    # center the covariates (or not)
    if(stdz==TRUE){
      X <- apply(X,2,function(x){(x-mean(x))/sd(x)})
      Z <- apply(Z,2,function(x){(x-mean(x))/sd(x)})
    }

    # fit the aft model
    if(corstr=="independence"){corstr1<-"ind"}else if(corstr=="exchangeable"){corstr1<-"ex"}
    res <- SMC.AFT.GEE(tobs=tobs,delta=delta,id=id,X=X,Z=Z,corstr=corstr1,
                       Var=Var,nboot=nboot,Intercept=TRUE,
                       nwrong.rate=5,control=SMC.AFT.GEE.control(
                         emmax=esmax,emeps=eps,emtrace=FALSE,
                         maxiter.latency=15,eps.latency = 1e-04,trace.latency = FALSE
                       )
    )[[2]]

  }

  ## tidy the results
  fit <- list()
  fit$call <- call
  fit$model <- model

  if(model=="ph"){

    if(Var==TRUE){

      # incidence part
      temp.incidence <- array(res$gamma, c(length(res$gamma), 4))
      rownames(temp.incidence) <- res$gamma_name
      colnames(temp.incidence) <- c("Estimate", "Std.Error",
                                    "Z value", "Pr(>|Z|)")
      temp.incidence[, 2] <- res$boots_gamma_sd
      temp.incidence[, 3] <- res$boots_gamma_zvalue
      temp.incidence[, 4] <- res$boots_gamma_pvalue
      # latency part
      temp.latency <- array(res$beta, c(length(res$beta), 4))
      rownames(temp.latency) <- res$beta_name
      colnames(temp.latency) <- c("Estimate", "Std.Error",
                                  "Z value", "Pr(>|Z|)")
      temp.latency[, 2] <- res$boots_beta_sd
      temp.latency[, 3] <- res$boots_beta_zvalue
      temp.latency[, 4] <- res$boots_beta_pvalue
      # for Estimated Correlation Parameters
      temp.cor <- array(0, c(2, 4))
      rownames(temp.cor) <- c("cor_incidence", "cor_latency")
      colnames(temp.cor) <- c("Estimate","Std.Error",
                              "Z value", "Pr(>|Z|)")
      temp.cor[, 1] <- c(res$alpha, res$rho)
      temp.cor[, 2] <- c(res$boots_gcor_sd, res$boots_bcor_sd)
      temp.cor[, 3] <- c(res$boots_gcor_zvalue, res$boots_bcor_zvalue)
      temp.cor[, 4] <- c(res$boots_gcor_pvalue, res$boots_bcor_pvalue)

    }else{

      # incidence part
      temp.incidence <- array(res$gamma, c(length(res$gamma), 1))
      rownames(temp.incidence) <- res$gamma_name
      colnames(temp.incidence) <- "Estimate"
      # latency part
      temp.latency <- array(res$beta, c(length(res$beta), 1))
      rownames(temp.latency) <- res$beta_name
      colnames(temp.latency) <- "Estimate"
      # for Estimated Correlation Parameters
      temp.cor <- array(0, c(2, 1))
      rownames(temp.cor) <- c("cor_incidence", "cor_latency")
      colnames(temp.cor) <- "Estimate"
      temp.cor[, 1] <- c(res$alpha, res$rho)

    }


  }else if(model=="aft"){

    if(Var==TRUE){

      # incidence part
      temp.incidence <- res$incidence
      colnames(temp.incidence) <- c("Estimate", "Std.Error",
                                    "Z value", "Pr(>|Z|)")
      temp.incidence[, 2] <- sqrt(res$incidence[,2])
      # latency part
      temp.latency <- res$latency
      colnames(temp.latency) <- c("Estimate", "Std.Error",
                                  "Z value", "Pr(>|Z|)")
      temp.latency[, 2] <- sqrt(res$latency[,2])
      # for Estimated Correlation Parameters
      temp.cor.est <- res$corr[,1]
      temp.cor.se  <- sqrt(res$corr[,2])
      temp.cor.zvalue <- temp.cor.est/temp.cor.se
      temp.cor.pvalue <- 2*(1-pnorm(abs(temp.cor.zvalue)))
      temp.cor <- data.frame(temp.cor.est,temp.cor.se,temp.cor.zvalue,temp.cor.pvalue)
      rownames(temp.cor) <- c("cor_incidence", "cor_latency")
      colnames(temp.cor) <- c("Estimate","Std.Error",
                              "Z value", "Pr(>|Z|)")

    }else{

      # incidence part
      temp.incidence <- res$incidence
      colnames(temp.incidence) <- "Estimate"
      # latency part
      temp.latency <- res$latency
      colnames(temp.latency) <- "Estimate"
      # for Estimated Correlation Parameters
      temp.cor <- res$corr[,1,drop=F]
      rownames(temp.cor) <- c("cor_incidence", "cor_latency")
      colnames(temp.cor) <- c("Estimate")

    }
    temp.latency <- temp.latency[-1,,drop=FALSE]

  }

  fit$incidence <- temp.incidence
  fit$latency   <- temp.latency
  fit$cor <- temp.cor
  fit$corstr <- corstr
  fit$num_of_clusters <- length(uid)
  fit$max_cluster_size <- max(table(id))

  class(fit) <- c("smgeecure")
  fit
  # print.smgeecure(fit)

}


#==== The main function for smc aft/ph mixture cure fitting model ====#
#' @title Print smgeecure object
#'
#' @description Output of \code{smgeecure} object.
#'
#' @aliases print.smgeecure
#'
#' @param fit an object of smgeecure
#'
#' @export print.smgeecure
print.smgeecure <- function(fit){

  ## print the fitting results
  if (!is.null(cl <- fit$call)) {
    cat("Call:\n")
    dput(cl)
  }

  cat("\nCure Probability Model:\n")
  print(fit$incidence)

  cat("\nFailure Time Distribution Model:\n")
  print(fit$latency)

  cat("\nEstimated Correlation Parameters:\n")
  print(fit$cor)

  cat("\nNumber of clusters:", fit$num_of_clusters)
  cat("\t Maximum cluster size:", fit$max_cluster_size)
  cat("\n")
  invisible(fit)

}


#========================================================#
# Fit the marginal semiparametric accelerated failure time mixture cure model ####
#========================================================#

#==== The main function for fitting model ====#
SMC.AFT.GEE <- function(tobs,delta,id,X,Z,corstr= c("ind","ar1","ex"),
                        Var=TRUE,nboot=100,nwrong.rate=0.4,Intercept=TRUE,
                        control = SMC.AFT.GEE.control()){

  #### Preparation: Basic Settings ####
  pX <- ncol(X); pZ <- ncol(Z)
  K <- length(unique(id)) # number of clusters
  uid <- sort(unique(id)) # unique ids (length(uid) is the number of individuals)

  ### Ests for original estimators (using all data)
  smcfit <- SMC.AFT.GEE.fit(tobs,delta,id,X,Z,corstr,Intercept=Intercept,control)

  ### bootstrap procedure (really time consuming!)
  the.boot.all <- array(NA,dim=c(nrow(smcfit$fit),ifelse(Intercept,pX+pZ+2+2,pX+pZ+1+2),nboot))
  if(Var==TRUE){
    iboot <- 1; nwrong <- 0
    while(iboot <= nboot){ # iboot <- 1
      cat("Current bootstrap number: ",iboot,"\n")
      repeat{
        if(nwrong > (nboot*nwrong.rate)){stop("An error!")}
        if("ex" %in% corstr | "ar1" %in% corstr){
          id.select.iboot <- sample((1:K), replace = TRUE)
          idx.iboot <- unlist(lapply(id.select.iboot,function(iid){which(id == iid)}))
          id.iboot <- unlist(lapply(1:K,function(ii){rep(ii,sum(id == id.select.iboot[ii]))}))
        }else{
          idx.iboot <- sample((1:length(tobs)), replace = TRUE)
          id.iboot <- 1:length(tobs)
        }
        ### fit the current models
        smctry <- try({
          smcfit.iboot <- SMC.AFT.GEE.fit(tobs[idx.iboot],delta[idx.iboot],id.iboot,
                                          X[idx.iboot,,drop=F],Z[idx.iboot,,drop=F],corstr,
                                          Intercept=Intercept,control)
        }, silent = T)
        if(is(smctry, "try-error") == TRUE | max(abs(smcfit.iboot$fit))>1e6){  #  | smcfit.iboot$convergence[1]==FALSE
          nwrong<-nwrong+1
          next
        }else{
          break
        }
      }
      ### Get the value for beta and gam
      the.boot.all[,,iboot] <- smcfit.iboot$fit
      ### next iteration
      iboot <- iboot + 1
    }
  }
  smcfit.var <- apply(the.boot.all,c(1,2),var)
  rm(the.boot.all)

  ### do inference and combine results
  if(Var==TRUE){mycol<-TRUE}else{mycol<-1}
  out0 <- rep(list(list()),nrow(smcfit$fit))
  for(imethod in 1:nrow(smcfit$fit)){
    # for latency part
    the.c <- smcfit$fit[imethod,]
    the.c.var <- smcfit.var[imethod,]
    zvalue.the.c <- the.c/sqrt(the.c.var)
    pvalue.the.c <- 2*(1-pnorm(abs(zvalue.the.c)))
    # combine
    idx.gam <- c(1:(pZ+1))
    if(Intercept==TRUE){
      idx.bet <- c((1:(pX+1))+(pZ+1))
      Intercept.Latency <- "(Intercept)"
    }else{
      idx.bet <- c((1:pX)+(pZ+1))
      Intercept.Latency <- NULL
    }

    out0[[imethod]] <- list(
      incidence = data.frame(Est=the.c[idx.gam],Var=the.c.var[idx.gam],zvalue=zvalue.the.c[idx.gam],
                             pvalue=pvalue.the.c[idx.gam],row.names=c('(Intercept)',colnames(Z)))[,mycol,drop=F],
      latency = data.frame(Est=the.c[idx.bet],Var=the.c.var[idx.bet],zvalue=zvalue.the.c[idx.bet],
                           pvalue=pvalue.the.c[idx.bet],row.names=c(Intercept.Latency,colnames(X)))[,mycol,drop=F],
      corr = data.frame(Est=the.c[-c(1:(pX+pZ+ifelse(Intercept,2,1)))],Var=the.c.var[-c(1:(pX+pZ+ifelse(Intercept,2,1)))])[,mycol,drop=F],
      convergence = smcfit$convergence[imethod]
    )
  }
  names(out0) <- c("fit.smcure",paste("fit.aftgee_",corstr,sep=""))

  ### extract output values
  out <- c(out0)
  return(out)

}


#==== The main function for fitting model ====#
SMC.AFT.GEE.fit <- function(tobs, delta, id, X, Z, corstr= c("ind","ar1","ex"),Intercept=FALSE,
                            control = SMC.AFT.GEE.control()){

  #### Preparation: Basic Settings ####
  pX <- ncol(X); pZ <- ncol(Z)
  nobs <- length(tobs)
  K <- length(unique(id)) # number of clusters
  uid <- sort(unique(id)) # unique ids (length(uid) is the number of individuals)
  n <- rep(0, length(uid)) # number of individuals in each cluster
  for (i in 1:length(uid)) {
    n[i] <- sum(id == i)
  }
  Xbar <- apply(X,2,mean)
  ZI <- as.matrix( as.data.frame(cbind(1, Z)) ) # with intercept 1
  XI <- as.matrix( as.data.frame(cbind(1, X)) ) # with intercept 1

  #### Give initial values before ES ####
  smcfit <- SMC.AFT.fit(tobs,delta,X,Z,em.maxit=control$emmax.smcure,em.eps=control$emeps)

  # dat.smcure <- data.frame(Time=tobs/100,delta,X)
  # pd <- smcure(Surv(Time,delta)~Z1+Z2+Z3+Z4+Z5,
  #              cureform=~Z1+Z2+Z3+Z4+Z5,data=dat.smcure,model="aft",
  #              Var = FALSE)
  # printsmcure(pd,Var = FALSE)
  #
  # data(e1684)
  # pd <- smcure(Surv(FAILTIME,FAILCENS)~TRT+SEX+AGE,
  #              cureform=~TRT+SEX+AGE,data=e1684,model="aft",
  #              Var = FALSE)
  # printsmcure(pd,Var = FALSE)

  if(Intercept==TRUE){
    beta.ori <- smcfit$latency[,1]
  }else{
    beta.ori <- smcfit$latency[-1,1]
  }
  gamma.ori <- smcfit$incidence[,1]
  fit.smcure <- c(gamma.ori,beta.ori)
  convergence.smcure <- smcfit$convergence

  #### Start my ES method #### (we need to convert it into a single function)
  fit.all <- array(0,dim=c(length(corstr), ifelse(Intercept,pX+pZ+2+2,pX+pZ+1+2) ))
  convergence.all <- rep(NA,length(corstr))
  for(icorstr in 1:length(corstr)){ # icorstr <- 1

    # assign initial values
    beta <- beta.ori
    gamma <- gamma.ori
    Shat <- SMC.AFT.St(tobs,delta,X,beta,delta,intercept=Intercept)$St
    iem <- 1
    repeat{

      ####################################
      ### update gij ( Y ) using old S ###
      ####################################
      expZgamma <- ZI %*% gamma
      Pic <- exp( expZgamma ) / ( 1 + exp( expZgamma ) )
      gm <- as.vector( delta + ( (1 - delta) * Pic * Shat ) /  ( 1 - Pic + Pic * Shat ) )

      ###########################################
      # obtain parameters in the incidence part #
      ###########################################
      incidence.fitm <- eval(parse(text = paste("geepack::geese", "(", "gm ~ Z", ", id = ", "id", ", family = binomial", ", corstr = '",
                                                corstr[icorstr], "'", ")", sep = "")))
      gamma.new <- incidence.fitm$beta
      # phi <- incidence.fitm$gamma
      # alpha <- incidence.fitm$alpha



      # latency.fitm <- eval(parse(text = paste("geepack::geese", "(", "yhat ~ X", ", id = ", "id", ", family = gaussian", ", corstr = '",
      #                                         corstr[icorstr], "'", ")", sep = "")))
      # beta.new <- latency.fitm$beta
      # beta.new <- latency.fitm$beta
      # # phi <- latency.fitm$gamma
      # # alpha <- latency.fitm$alpha


      #########################################
      # obtain parameters in the latency part #
      #########################################
      ### initial value for this iteration using Gehan's method (Zhang and Peng 2007)
      beta0 <- optim(par=rep(0,ifelse(Intercept,pX+1,pX)),fn=SMC.AFT.rank,method="Nelder-Mead",
                     control=list(maxit=500,reltol=1e-04),
                     yobs=tobs,delta=delta,X=X,w=gm,intercept=Intercept)$par
      # as initial beta has been given: beta0 (we will use GEE to make it better)
      for(iupdate in 1:control$maxiter.latency){
        # previous beta
        betaprev <- beta0

        # some values
        eres <- SMC.AFT.eRes(tobs,delta,X,beta0,gm,intercept=Intercept)$eres
        if(Intercept==TRUE){
          yhat <- delta * log(tobs) + (1 - delta) * (eres + XI %*% beta0)
        }else{
          yhat <- delta * log(tobs) + (1 - delta) * (eres + X %*% beta0)
        }


        # # repeat until rho converges
        # beta.old <- beta0
        # repeat{

        # calculate current Omega, phi and rho
        if(Intercept==TRUE){
          r.hat <- as.vector(yhat - XI %*% beta0)
        }else{
          r.hat <- as.vector(yhat - X %*% beta0)
        }
        phi.hat <- sum( r.hat^2 ) / ( nobs - ifelse(Intercept,pX+1,pX) )
        # calculate rho
        r.m <- matrix(0, nrow = K, ncol = max(n)) # change the format of r.hat (K*max(n))
        for (i in 1:K) {
          r.m[i,1:n[i]] <- r.hat[id == i]
        }
        if(corstr[icorstr]=="ind"){
          # for ind
          rho <- 0
        }else if(corstr[icorstr]=="ex"){
          # for ex
          rho0 <- 0
          for (i in 1:K) {
            if (n[i] == 1) {
              rho0 <- rho0 + r.m[i, 1]
            } else {
              for (j in 1:(n[i] - 1)) rho0 <- rho0 + r.m[i, j] * sum(r.m[i,(j + 1):n[i]])
            }
          }
          rho <- (phi.hat^(-1)) * rho0/(sum(n * (n - 1))/2 - ifelse(Intercept,pX+1,pX) )
        }else if(corstr[icorstr]=="ar1"){
          # for ar1
          rho0 <- dividenumber <- rep(0,max(n))
          for(i in 1:K){
            for(k in 0:(n[i]-1)){
              for(j in 1:(n[i]-k)){
                rho0[k+1] <- rho0[k+1] + r.m[i, j] * r.m[i, j+k]
                dividenumber[k+1] <- dividenumber[k+1] + 1
              }
            }
          }
          dividenumber <- dividenumber - ifelse(Intercept,pX+1,pX)
          rho <- (phi.hat^(-1)) * rho0 / dividenumber
          # rho[3] <- rho[2]^2
        }

        # solve (for the current rho)
        I1 <- array(0,dim=c(ifelse(Intercept,pX+1,pX),ifelse(Intercept,pX+1,pX)))
        I2 <- array(0,dim=c(ifelse(Intercept,pX+1,pX),1))
        for (j in 1:K) {
          Omega <- CWM.corstr(rho,n[j],corstr[icorstr]) # choose appropriate type of Omega
          Omega.inv <- solve(Omega*phi.hat)
          K.yhatj <- yhat[id==j]
          if(Intercept==TRUE){
            mu.Xj <- XI[id==j,,drop=FALSE]
            gmj <- gm[id==j]
            I1j <- t(mu.Xj) %*% Omega.inv %*% diag(gmj,nrow=n[j],ncol=n[j]) %*% mu.Xj
            I2j <- t(mu.Xj) %*% Omega.inv %*% diag(gmj,nrow=n[j],ncol=n[j]) %*% K.yhatj
          }else{
            mu.Xj <- X[id==j,,drop=FALSE]
            gmj <- gm[id==j]
            I1j <- (t(mu.Xj)-Xbar) %*% Omega.inv %*% diag(gmj,nrow=n[j],ncol=n[j]) %*% mu.Xj
            I2j <- (t(mu.Xj)-Xbar) %*% Omega.inv %*% diag(gmj,nrow=n[j],ncol=n[j]) %*% K.yhatj
          }

          I1 <- I1 + I1j
          I2 <- I2 + I2j
        }
        beta0 <- as.vector( solve(I1) %*% I2 )

        #   # whether to stop
        #   break
        #   if(max(abs(beta.old-beta0))>1e-6){
        #     beta.old <- beta0
        #   }else{
        #     break
        #   }
        #
        # }

        if (control$trace.latency) cat("\n beta:", as.numeric(beta0), "\n")
        if (max(abs(beta0 - betaprev) / abs(beta0)) <= control$eps.latency) break

      } # => get a new beta (beta)
      beta.new <- beta0

      ### update Shat
      Shat.new <- SMC.AFT.St(tobs,delta,X,beta.new,gm,intercept=Intercept)$St

      ########################
      # when to break our ES #
      ########################
      convergence.value <- sum(c(gamma.new - gamma, beta.new - beta)^2)
      if (control$emtrace){
        print('-----')
        print(iem)
        print(drop(c(gamma.new, beta.new)))
        print(convergence.value)
      }
      if ( ( convergence.value > control$emeps ) && ( iem < control$emmax ) ) {
        gamma <- gamma.new
        beta <- beta.new
        Shat <- Shat.new
        iem <- iem + 1
      } else break
    }# end all ES
    convergence <- (iem<control$emmax & convergence.value < control$emeps)

    # assgin fitted values
    if(corstr[icorstr]=='ind'){
      alp <- rho <- 0
    }else{
      alp <- incidence.fitm$alpha
      if(corstr[icorstr]=='ar1'){rho <- rho[2]}
    }
    fit.all[icorstr,] <- c(gamma.new,beta.new,alp,rho)
    convergence.all[icorstr] <- convergence

  } # end corstr


  ######################
  # Output our results #
  ######################
  out <- list(fit=rbind(c(fit.smcure,0,0),fit.all),
              convergence=c(convergence.smcure,convergence.all))
  Intercept.Latency <- ifelse(Intercept,"Intercept",NULL)
  colnames(out$fit) <- c("Intercept",paste("gamma",1:pZ, sep=""),
                         Intercept.Latency,paste("beta",1:pX, sep=""),
                         "cor_incidence","cor_lantency")
  rownames(out$fit) <- c("smcure",paste("aftgee_",corstr,sep=""))
  out

  return(out)

}


#====  Function used to control some basic settings ====#
SMC.AFT.GEE.control <- function(emmax=15,
                                emeps=1e-04,
                                emtrace=FALSE,
                                emmax.smcure=50,
                                maxiter.latency = 15,
                                eps.latency = 1e-04,
                                trace.latency = FALSE) {
  list(emmax=emmax, emeps=emeps,emtrace = emtrace,
       emmax.smcure=emmax.smcure,
       maxiter.latency = maxiter.latency, eps.latency = eps.latency, trace.latency = trace.latency)

}





#========================================================#
# Fit the semiparametric accelerated failure time mixture cure model using Zhang and Peng (2007) ####
#========================================================#

#==== The main function for smc aft mixture cure fitting model ====#
SMC.AFT.fit <- function(yobs,delta,X,Z,incidence=c("logit"),em.maxit=50,em.eps=1e-4){

  ### Specify the dimension and prepare data
  N <- length(yobs) # number of individuals
  pX <- ncol(X)+1
  pZ <- ncol(Z)+1
  ZI <- as.matrix(cbind(rep(1,N),Z)) # [dataframe: for incidence part]

  ### do EM algorithm
  # obtain initial bet, gam and baseline nonparametric part
  bet.old <- survival::survreg(survival::Surv(yobs,delta)~X)$coef # from survival package
  gam.old <- eval(parse(text=paste("glm", "(", "delta~Z",",family = quasibinomial(link='",incidence,"'",")",")",sep = "")))$coef
  St.old <- SMC.AFT.St(yobs,delta,X,bet.old,delta)$St
  numit <- 1
  repeat{

    # calculate w first
    expZgam <- exp(ZI%*%gam.old)
    uncureprob <- expZgam/(1+expZgam)
    w <- delta+(1-delta)*(uncureprob*St.old)/(1-uncureprob+uncureprob*St.old)

    # update gam
    gam <- eval(parse(text=paste("glm","(","w~Z",",family = binomial(link='",incidence,"'", ")",")", sep = "")))$coef
    # update bet
    bet <- optim(par=rep(0,pX),fn=SMC.AFT.rank,# method="BFGS",
                 # control=list(maxit=500,reltol=em.eps),
                 yobs=yobs,delta=delta,X=X,w=w)$par
    # update St
    St <- SMC.AFT.St(yobs,delta,X,bet,w)$St

    ### update or stop
    convergence.value <- sum(c(gam-gam.old,bet-bet.old)^2) #+ sum((St-St.old)^2)
    if(convergence.value >= em.eps & numit < em.maxit){
      bet.old <- bet; gam.old <- gam
      St.old <- St
      numit <- numit + 1
    }else{
      break
    }

  } # end em reps
  convergence <- (numit<em.maxit & convergence.value < em.eps)

  ### extract values
  out <- list(
    latency = data.frame(Est=bet,row.names=c('Intercept',colnames(X))),
    incidence = data.frame(Est=gam,row.names=c('Intercept',colnames(Z))),
    convergence = convergence
  )
  return(out)

}


#==== the objective funtion for rank-based method ====#
SMC.AFT.rank <- function(bet,yobs,delta,X,w,intercept=TRUE){

  if(intercept==FALSE){Xbet<-as.vector(X%*%bet)}else{Xbet<-as.vector(cbind(1,X)%*%bet)}
  error <- log(yobs) - Xbet
  LossGehan <- mean(sapply(1:length(yobs),function(i){
    sum((error[i]<error)*abs(error-error[i])*w*delta[i])}))
  return(LossGehan)

}



#==== estimate my S(t|X) in each step during em ====#
SMC.AFT.St <- function(yobs,delta,X,bet,w,intercept=TRUE){
  # estimate my S(t|X) in each step during em

  # prepare error
  if(intercept==FALSE){Xbet<-as.vector(X%*%bet)}else{Xbet<-as.vector(cbind(1,X)%*%bet)}
  error <- log(yobs) - Xbet
  # calculate Hazards
  sumRisk <- sapply(error,function(errori){sum((error>=errori)*w)})
  ht <- ifelse(delta==1,delta/sumRisk,0)
  Ht <- sapply(error,function(tmi){sum((error<=tmi)*ht)})
  Ht[error>max(error[delta==1])] <- Inf
  Ht[error<min(error[delta==1])] <- 0
  # calculate Survivals
  out <- list(St = exp(-Ht))
  # output
  return(out)
}


#====  Function used to Calculate Yhat,shat in each step during em ====#
SMC.AFT.eRes <- function(tobs,delta,X,beta,gm=rep(1,length(yobs)),intercept=TRUE) {
  # Function used to Calculate Yhat in each step during em

  if(intercept==FALSE){Xbet<-as.vector(X%*%beta)}else{Xbet<-as.vector(cbind(1,X)%*%beta)}
  error <- log(tobs) - Xbet
  nobs <- length(error)
  ord <- order(error)
  errori <- error[ord]
  deltai <- delta[ord]
  gmi <- gm[ord]
  Shati <- SMC.AFT.St(tobs,delta,X,beta,gm,intercept)$St[ord]
  # Shati <- approx(errori, Shati, errori)$'y'

  # calculate 1
  errordif <- c(diff(errori), 0)   ## diff(ei) gives 1 less terms
  ehat <- rev(cumsum(rev(errordif * Shati)))  # round(ehat,5)
  ehat <- ehat/Shati + errori    ## + errori because there was a diff() in errordif
  ehat[is.na(ehat)] <- 0#errori[is.na(ehat)]
  eres <- ehat
  eres[c(1:nobs)[ord]] <- ehat  ## puting it back to the original order
  # # # calculate 2
  # # ehat2 <- rev(cumsum(rev(errori * errordif * Shati)))
  # # ehat2 <- 2 * ehat2/Shati + errori^2
  # # ehat2[is.na(ehat2)] <- errori[is.na(ehat2)]^2
  # # ehat2[which(ehat2 < 0)] <- NaN
  # # eres2 <- ehat2
  # # eres2[c(1:nobs)[ord]] <- ehat2

  # Shatidif <- -c(diff(c(Shati)),0)   ## diff(ei) gives 1 less terms
  # ehat <- rev(cumsum(rev(errori * Shatidif)))  # round(ehat,5)
  # ehat <- ehat/Shati
  # ehat[is.na(ehat)] <- 0
  # eres <- ehat
  # eres[c(1:nobs)[ord]] <- ehat  ## puting it back to the original order


  # output
  out <- list(eres=eres)
  return(out)
}




#====  Function used to construct the correlation working matrix ====#
CWM.corstr <- function(rho,ndim,corstr){
  # construct correlation working matrix

  if(corstr=='ind'){
    CWM <- matrix(0, ndim, ndim)
  }else if(corstr=='ex'){
    CWM <- matrix(rho, ndim, ndim)
  }else if(corstr=='ar1'){
    CWM <- matrix(0,ndim,ndim)
    # for(i in 1:ndim){for(j in i:ndim){CWM[i,j] <- CWM[j,i] <- rho^{abs(i-j)}}}
    # for(i in 1:ndim){for(j in i:ndim){CWM[i,j] <- CWM[j,i] <- rho[j-i+1]}}
    for(i in 1:ndim){for(j in i:ndim){CWM[i,j] <- CWM[j,i] <- rho[2]^(abs(i-j))}}
  }
  diag(CWM) <- 1

  return(CWM)
}



#========================================================#
# The following codes contain all the documentation of datasets  ####
#========================================================#


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




