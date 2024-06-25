# smgeecure
This is an R package for fitting marginal semiparametric mixture cure models for clustered survival data.
- Two semiparametric mixture cure models are implemented: the *marginal accelerated failure time mixture cure* (**AFTMC**) model and the *marginal proportional hazards mixture cure* (**PHMC**) model.


## Package description and included main functions

Our provided main functions are (we refer to their help pages for more details):
- **smgeecure**: fit the models in various ways with syntax
```R
SMC.AuxSP(formula, cureform, sdata, aux = NULL, hetero = FALSE, N = Inf, latency = "PH", nboot = 400)
```
- **print.smgeecure**: print outputted results from SMC.AuxSP() with syntax
```R
print.SMC.AuxSP(object)
```

## Two numerical illustrations

### An example using a simulated dataset is shown below:


- Two real data examples are presented in the help file of *smgeecure()*.

