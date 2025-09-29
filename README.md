# NCVMP for MELS

```
library(lme4)
source("mels_ncvmp.R")
source("fast_regls.R")

##################################
### Example 1. Riesby Data Set ###
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
