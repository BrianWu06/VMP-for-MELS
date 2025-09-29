# NCVMP for MELS

This repo is an implementation of a Non-Conjugate Variational Message Passing algorithm for the Mixed-Effect Location-Scale Model. 

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

The summary should look lile: 
```
> summary(riesby_ncvmp_r)
## NCVMP for MELS ##
--------------------------------------------------------
--- Mean Model Parameters (beta) ---
            Estimate Std. Error  z value p-value
(Intercept)  22.3375     0.3903  57.2268  0.0000
week         -2.3072     0.1522 -15.1549  0.0000
endog         1.8771     0.5681   3.3043  0.0010
endweek      -0.0276     0.2231  -0.1239  0.9014

--- Between-Subject Variance Parameters (alpha) ---
            Estimate Std. Error z value p-value
(Intercept)   2.2449     0.2627   8.547  0.0000
endog         0.5069     0.3508   1.445  0.1485

--- Within-Subject Variance Parameters (tau) ---
            Estimate Std. Error z value p-value
(Intercept)   2.1324     0.1481 14.3973  0.0000
week          0.1916     0.0413  4.6420  0.0000
endog         0.2823     0.1467  1.9244  0.0543

--- Random Effect Standard Deviation ---
  Std. Dev of omega (scale): 0.5757

--- Mean/Variance Interaction Parameter (r) ---
  Estimate Std. Error z value p-value
r   0.0453     0.0207  2.1892  0.0286
-------------------------------------------------------
Convergence Details:
  Algorithm converged in 24 iterations.
  Total Runtime: 0.11 secs
```
