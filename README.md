# smgeecure
This is an R package for fitting marginal semiparametric mixture cure models for clustered survival data.
- Two semiparametric mixture cure models are implemented: the *marginal accelerated failure time mixture cure* (**AFTMC**) model and the *marginal proportional hazards mixture cure* (**PHMC**) model.

## Package description and included main functions

Installation of this package can be done locally after downloading the package manually from this github website. We have also uploaded this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package using R command
```R
install.packages("smgeecure").
```

Our provided main functions are (we refer to their help pages for more details):
- **smgeecure**: fit the models in various ways with synopsis
```R
smgeecure(formula, cureform, id, data, model = c("aft", "ph"), corstr = c("independence", "exchangeable", "ar1"),
          Var = TRUE,nboot = 100, stdz = TRUE, esmax = 20, eps = 1e-04)
```
- **print.smgeecure**: print outputted results from **smgeecure()** with syntax
```R
print.SMC.AuxSP(object)
```

## Numerical illustrations

### An example using a real dataset from TCGA program is shown below:
```R
# ---- prepare the data
data(tonsil)
tonsil <- tonsil[-c(141,136,159),]
tonsil$Sex <- ifelse(tonsil$Sex == 1, 0, 1) # 1="Female"
tonsil$Cond <- ifelse(tonsil$Cond == 1, 0, 1) # 0=no disability
tonsil$T <- ifelse(tonsil$T < 4, 0, 1)
tonsil$Grade2 <- ifelse(tonsil$Grade==2,1,0)
tonsil$Grade3 <- ifelse(tonsil$Grade==3,1,0)
table(tonsil$Inst)

# ---- plot the KM survival curve (overall)
pdf("Figure_Tonsil_KM.pdf",width=7,height=6)
plot(survival::survfit(survival::Surv(Time, Status) ~ 1, data = tonsil), 
     conf.int = T, mark.time = TRUE, lwd=1.5,
     ylab = "Survival Probability", xlab = "Survival Time (in Days)", 
     xlim = c(0,2000), ylim = c(0,1),
     cex.lab = 1, 
     xaxt = "n", yaxt = "n"
)
axis(1,seq(0,2000,500),seq(0,2000,500),cex.axis = 1)
axis(2,seq(0,1,0.2),seq(0,1,0.2),cex.axis = 1)
dev.off()

# ---- plot the KM survival curve (stratified)
pdf("Figure_Tonsil_KM_SexTumorsize.pdf",width=7,height=6)
plot(survival::survfit(survival::Surv(Time, Status) ~ Sex + T, data = tonsil),
     conf.int = F,  mark.time = TRUE, lwd=1.5,
     ylab = "Survival Probability", xlab = "Survival Time (in Days)",
     xlim = c(0,2000), ylim = c(0,1),
     lty = c(1,2,3,4), col = c(1,2,3,4),
     xaxt = "n", yaxt = "n", cex.lab=1)
axis(1,seq(0,2000,400),seq(0,2000,400),cex.axis=1)
axis(2,seq(0,1,0.2),seq(0,1,0.2),cex.axis=1)
legend("topright", c("Massive/Female","Massive/Male","Non-massive/Female","Non-massive/Male"),
       cex = 1, col = c(1,2,3,4), lty = c(1,2,3,4))
dev.off()

# ---- plot the log-cumulative survival curve (check ph assumption)
#         Plot BY Sex 
# estimates for 1
KM.fit.1 <- survival::survfit(survival::Surv(Time, Status) ~ 1,
                               data = tonsil, subset = (Sex==1))
Infty.1 <- min(KM.fit.1$surv)
surv.1.Latency <- (KM.fit.1$surv-Infty.1) / (1-Infty.1)
time.1.Latency <- KM.fit.1$time
cumu.1.Latency.Nozero <- -log(surv.1.Latency[surv.1.Latency!=0])
time.1.Latency.Nozero <- KM.fit.1$time[surv.1.Latency!=0]
# estimates for 0
KM.fit.0 <- survival::survfit(survival::Surv(Time, Status) ~ 1,
                               data = tonsil, subset = (Sex==0))
Infty.0 <- min(KM.fit.0$surv)
surv.0.Latency <- (KM.fit.0$surv-Infty.0) / (1-Infty.0)
time.0.Latency <- KM.fit.0$time
cumu.0.Latency.Nozero <- -log(surv.0.Latency[surv.0.Latency!=0])
time.0.Latency.Nozero <- KM.fit.0$time[surv.0.Latency!=0]
# plot of log-cumulative function
pdf("Figure_Tonsil_LogCum_Sex.pdf",width=7,height=6)
plot(log(cumu.1.Latency.Nozero)~log(time.1.Latency.Nozero),type="s",lty=1,col="red",
     xlim = c(0,8), ylim=c(-5,1.5), lwd=1.5,
     ylab = "Logarithm of the Cumulative Hazard Function", 
     xlab = "Logarithm of the Survival Time",
     xaxt = "n",
     cex.lab=1
)
lines(log(cumu.0.Latency.Nozero)~log(time.0.Latency.Nozero),type="s",lty=2,col="black",lwd=1.5)
axis(1,seq(0,8,2),seq(0,8,2),cex.axis = 1)
legend("bottomright", c("Female","Male"), cex = 1, 
       col = c("red","black"), lty = c(1,2))
dev.off()

# ---- plot the log-cumulative survival curve (check ph assumption)
#         Plot BY Sex 
#         (Just use uncensored patients)
UncuredPatients <- tonsil$Status==1 
pdf("Figure_Tonsil_LogCum_Uncen.pdf",width=7,height=6)
xlim_lower <- 0; xlim_upper <- 1800; xlim_step  <- 300
ylim_lower <- -5; ylim_upper <- 2; ylim_step  <- 1
mylty <- c(1,2); mycol <- c("black","red")
plot(survival::survfit(survival::Surv(Time, Status) ~ Sex, 
                       data = tonsil, subset=UncuredPatients,
                       stype=1), 
     conf.int = F,  mark.time = FALSE, cumhaz = TRUE, log="xy",
     ylab = "Logarithm of the Cumulative Hazard Function", 
     xlab = "Logarithm of the Survival Time (in Days)", 
     # xlim = c(xlim_lower,xlim_upper), ylim = c(ylim_lower,ylim_upper),
     lty = mylty, col = mycol, lwd=1.5,
     xaxt = "n", yaxt = "n")
axis(1,c(100,300,600,1500),c(100,300,600,1500))
axis(2,exp(seq(ylim_lower,ylim_upper,ylim_step)),seq(ylim_lower,ylim_upper,ylim_step))
legend("bottomright", c("male","female"), cex = 1.3, #bty="n",
       col = mycol, lty = mylty)
dev.off()

# ---- Fit the data using marginal semi-parametric marginal AFTMC model

# fit marginal semi-parametric marginal AFTMC model (exchangeable correlation)
set.seed(521)
tonsil.aft.gee.ex <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "exchangeable", Var = T, nboot = 100
)
print.smgeecure(tonsil.aft.gee.ex)

# fit marginal semi-parametric AFTMC model (independence correlation)
set.seed(521)
tonsil.aft.gee.ind <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "independence", Var = T, nboot = 100
)
print.smgeecure(tonsil.aft.gee.ind)

# output the results
write.csv(cbind(tonsil.aft.gee.ex$incidence,rbind(NA,tonsil.aft.gee.ex$latency)),
          "tonsil.aft.gee.ex.csv")
write.csv(cbind(tonsil.aft.gee.ind$incidence,rbind(NA,tonsil.aft.gee.ind$latency)),
          "tonsil.aft.gee.ind.csv")

# ---- Fit the data using marginal semi-parametric marginal PHMC model

# fit marginal semi-parametric PHMC model (exchangeable correlation)
set.seed(1)
tonsil.ph.gee.ex <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "ph", corstr = "exchangeable", Var = T,nboot = 100
)
print.smgeecure(tonsil.ph.gee.ex)

# fit marginal semi-parametric PHMC model (independence correlation)
set.seed(1)
tonsil.ph.gee.ind <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "ph", corstr = "independence", Var = F, nboot = 100
)
print.smgeecure(tonsil.ph.gee.ind)
```

More practical examples are presented in the help file of *smgeecure()*.

