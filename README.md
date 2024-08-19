# *smgeecure* R package
## Fitting Semiparametric Marginal Mixture Cure Models
This is an R package for fitting marginal semiparametric mixture cure models for clustered survival data.
- Two semiparametric mixture cure models are implemented: the *marginal accelerated failure time mixture cure* (**AFTMC**) model and the *marginal proportional hazards mixture cure* (**PHMC**) model.

## Package description and included main functions

Installation of this package can be done locally after downloading the package manually from this github website. We will also upload this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package. Currently, it can be loaded using R command
```R
devtools::install_github("yiniu06/smgeecure")
library(smgeecure)
```

The main function included in our R package is *smgeecure()* and there is also a function *print.smgeecure()* for printing fitted results with a better presentation. To sum up, they can be called via:
- **smgeecure**: fit the models in various ways with synopsis
```R
smgeecure(formula, cureform, data, id, model = c("aft", "ph"),
          corstr = c("independence", "exchangeable", "ar1"),
          Var = TRUE, nboot = 100, stdz = TRUE, esmax = 20, eps = 1e-04)
```
- **print.smgeecure**: print outputted results from the previous function *smgeecure()* with syntax
```R
print.smgeecure(fit)
```
We refer to their help pages for more detailed explanations of the corresponding arguments. Visually, the usages of part of the key arguments can be summarized via the following Figure:
![Graph-Illustration-of-Our-Method](https://github.com/user-attachments/assets/744d9d80-52af-459e-ac08-496be9486d0d)



## Numerical illustrations

An example of the use of this package using real data can be found in the following sections of this file.

### An example using the tonsil cancer clinical trial data

The tonsil cancer clinical trial study was conducted by the Radiation Therapy Oncology Group in the United States. We consider the semparametric marginal AFTMC model to the data with our function $smgeecure()$.

#### Data preparation
```R
data(tonsil)
tonsil <- tonsil[-c(141,136,159),]
tonsil$Sex <- ifelse(tonsil$Sex == 1, 0, 1) # 1="Female"
tonsil$Cond <- ifelse(tonsil$Cond == 1, 0, 1) # 0=no disability
tonsil$T <- ifelse(tonsil$T < 4, 0, 1)
tonsil$Grade2 <- ifelse(tonsil$Grade==2,1,0)
tonsil$Grade3 <- ifelse(tonsil$Grade==3,1,0)
table(tonsil$Inst)
```

#### Plot the stratified Kaplan-Meier survival curve
```R
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
```
![Figure_Tonsil_KM_SexTumorsize](https://github.com/user-attachments/assets/7874b5c8-46ea-4235-af6b-6e6e49c592ac)


#### Fit the data using marginal semi-parametric marginal AFTMC model
- exchangeable correlation
```R
set.seed(521)
tonsil.aft.gee.ex <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "exchangeable", Var = T, nboot = 100
)
print.smgeecure(tonsil.aft.gee.ex)
```
- independence correlation
```R
set.seed(521)
tonsil.aft.gee.ind <- smgeecure(
        formula = Surv(Time, Status) ~ Sex + factor(Grade) + Age + Cond + T, 
        cureform = ~ Sex + factor(Grade) + Age + Cond + T, id = tonsil$Inst, 
        data = tonsil, model = "aft", corstr = "independence", Var = T, nboot = 100
)
print.smgeecure(tonsil.aft.gee.ind)
```

### An example using the bone marrow transplantation data

The bone marrow transplantation data consists of 137 patients enrolled in the study from March 1, 1984 to June 30, 1989 at four institutions or hospitals. We consider the semparametric marginal PHMC model to the data with our function $smgeecure()$.

#### Data preparation
```R
data(bmt)
bmt$g <- factor(bmt$g, label = c("ALL", "AML low risk","AML high risk"))
bmt$Z8 <- factor(bmt$Z8, label = c("Otherwise", "FAB"))
```

#### Plot the stratified Kaplan-Meier survival curve
```R
plot(survival::survfit(survival::Surv(T2, d3) ~ factor(g), data = bmt),
     conf.int = F,  mark.time = TRUE, lwd=1.5,
     ylab = "Survival Probability", xlab = "Survival Time (in Days)",
     xlim = c(0,2800), ylim = c(0,1),
     lty = c(1,2,3), col = c(1,2,3),
     xaxt = "n", yaxt = "n", cex.lab=1)
axis(1,seq(0,2800,400),seq(0,2800,400),cex.axis=1)
axis(2,seq(0,1,0.2),seq(0,1,0.2),cex.axis=1)
legend("topright", c("ALL","AML low risk","AML high risk"),
       cex = 1, col = c(1,2,3), lty = c(1,2,3))
```
![Figure_bmt_KM_g](https://github.com/user-attachments/assets/e7a48faf-2e77-49aa-84cf-3190e15e0cae)


#### Fit the data using marginal semi-parametric marginal AFTMC model
- exchangeable correlation
```R
set.seed(1)
bmt.ph.gee.ex <- smgeecure(
  formula = Surv(T2, d3) ~ factor(g) + Z8, cureform = ~ factor(g) + Z8, 
  id = bmt$Z9, data = bmt, model = "ph", corstr = "exchangeable", 
  Var = T,nboot = 100, esmax = 100, eps = 1e-06
)
print.smgeecure(bmt.ph.gee.ex)
```
- independence correlation
```R
set.seed(1)
bmt.ph.gee.ind <- smgeecure(
  formula = Surv(T2, d3) ~ factor(g) + Z8, cureform = ~ factor(g) + Z8, 
  id = bmt$Z9, data = bmt, model = "ph", corstr = "independence", 
  Var = T, nboot = 100, esmax = 100, eps = 1e-06
)
print.smgeecure(bmt.ph.gee.ind)
```



