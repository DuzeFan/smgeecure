# smgeecure
This is an R package for fitting marginal semiparametric mixture cure models for clustered survival data.
- Two semiparametric mixture cure models are implemented: the *marginal accelerated failure time mixture cure* (**AFTMC**) model and the *marginal proportional hazards mixture cure* (**PHMC**) model.

## Package description and included main functions

Installation of this package can be done locally after downloading the package manually from this github website. We have also uploaded this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package using R command
```R
install.packages("smgeecure").
```

Our provided main functions are (we refer to their help pages for more details):
- **smgeecure**: fit the models in various ways with syntax
```R
smgeecure(formula, cureform, id, data, model = c("aft", "ph"), corstr = c("independence", "exchangeable", "ar1"),
          Var = TRUE,nboot = 100, stdz = TRUE, esmax = 20, eps = 1e-04)
```
- **print.smgeecure**: print outputted results from SMC.AuxSP() with syntax
```R
print.SMC.AuxSP(object)
```

## Numerical illustrations

### An example using a real dataset from TCGA program is shown below:
```R
library(smgeecure)

```

More practical examples are presented in the help file of *smgeecure()*.

