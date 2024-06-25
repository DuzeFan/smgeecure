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

