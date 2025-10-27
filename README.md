# VMP for MELS

This repo is an implementation of a Variational Message Passing algorithm for the Mixed-Effect Location-Scale Model. 

The MELS model is as follow:    
$y_{ij} = X_{ij}^\top \beta + \nu_i + \epsilon_{ij}$    
$\nu_i = \mathcal{N}(0, \exp{(u_i^\top \alpha)})$    
$\epsilon_{ij} = \mathcal{N}(0, \exp{(w_{ij}^\top \tau + \omega_i)})$    
$\omega_i = \mathcal{N}(0, \sigma^2_{\omega})$    

The mels_ncvmp takes the following input variables: 
1. y: a string that indicates the name of the response variable.
2. beta_formula: the mean model formula.
3. alpha_formula: the BS variance model formula, this can only allow for non-time varying variable for precise parameter estimation.
4. tau_formula: the WS variance model formula.
5. id: a string of the id variable's name.
6. data: a dataframe of the data set to be processed.
7. max_iter: the maximum number of iterations.
8. tol: the convergence criterion.
9. beta_prior, alpha_prior, tau_prior, A_prior: the priors for beta, alpha, tau and A espectively.
10. verbose: whether to print additional information or not.

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
                           id = "id"
                           data = riesby)
summary(riesby_ncvmp)

riesby_ncvmp_boot <- bootstrap_mels_ncvmp(riesby_ncvmp, B = 1000)
summary(riesby_ncvmp_boot)
```

The summary should look lile: 
```
> riesby_ncvmp <- mels_ncvmp(y = "hamd", 
+                            beta_formula = ~ week + endog + endweek, 
+                            alpha_formula = ~ endog, 
+                            tau_formula = ~ week + endog, 
+                            id = "id",
+                            data = riesby)
CAVI converges at iteration  55 
> summary(riesby_ncvmp)
## NCVMP for MELS ##
--------------------------------------------------------
--- Mean Model Parameters (beta) ---
            Estimate Std. Error  z value p-value
(Intercept)  22.2573     0.3884  57.3107  0.0000
week         -2.2673     0.1505 -15.0615  0.0000
endog         1.8679     0.5694   3.2808  0.0010
endweek      -0.0139     0.2221  -0.0626  0.9501

--- Between-Subject Variance Parameters (alpha) ---
            Estimate Std. Error z value p-value
(Intercept)   2.2777     0.2626  8.6732  0.0000
endog         0.4937     0.3507  1.4077  0.1592

--- Within-Subject Variance Parameters (tau) ---
            Estimate Std. Error z value p-value
(Intercept)   2.1381     0.1475 14.4928   0.000
week          0.1841     0.0411  4.4830   0.000
endog         0.2928     0.1467  1.9957   0.046

--- Random Effect Standard Deviation ---
  Std. Dev of omega (scale): 0.6101
-------------------------------------------------------
Convergence Details:
  Algorithm converged in 55 iterations.
  Total Runtime: 0.13 secs 
```
