# NCVMP for MELS

The MELS model is as follow: 
$y_{ij} = X_{ij}^\top \beta + \nu_i + \epsilon_{ij}$
$\nu_i = \mathcal{N}(0, \exp{u_i^\top \alpha})$
$\epsilon_{ij} = \mathcal{N}(0, \exp{w_{ij}^\top \tau} + \omega_i)$
$\omega_i = \mathcal{N}(0, \sigma^2_{\omega})$

riesby_ncvmp takes the following inputs: 
1. y: a string that indicates the name of the response variable.
2. beta_formula: the mean model formula.
3. alpha_formula: the BS variance model formula.
4. tau_formula: the WS variance model formula.
5. id: a string of the id variable's name.
6. r: whether to estimate the parameter r.
7. data: a dataframe of the data set to be processed.
8. max_iter: the maximum number of iterations.
9. tol: the convergence criterion.
10. beta_prior, alpha_prior, tau_prior, A_prior, r_prior: the priors for beta, alpha, tau, A and r respectively.
11. verbose: whether to print additional information or not.

```
library(lme4)
source("mels_ncvmp.R")
source("fast_regls.R")

##################################
### Example: Riesby Data Set #####
##################################

riesby <- read.table("RIESBY.DAT.txt", na.strings = ".")
colnames(riesby) <- c("id", "hamd", "intcpt", "week", "endog", "endweek")

riesby_ncvmp <- mels_ncvmp(y = "hamd", 
                           beta_formula = ~ week + endog + endweek, 
                           alpha_formula = ~ endog, 
                           tau_formula = ~ week + endog, 
                           id = "id", 
                           r = FALSE, 
                           data = riesby)
summary(riesby_ncvmp)

riesby_ncvmp_r <- mels_ncvmp(y = "hamd", 
                             beta_formula = ~ week + endog + endweek, 
                             alpha_formula = ~ endog, 
                             tau_formula = ~ week + endog, 
                             id = "id", 
                             r = TRUE, 
                             data = riesby)
summary(riesby_ncvmp_r)

riesby_ncvmp_boot <- bootstrap_mels_ncvmp(riesby_ncvmp, B = 1000)
summary(riesby_ncvmp_boot)
```
