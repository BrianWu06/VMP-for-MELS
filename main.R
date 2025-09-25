library(lme4)
source("mels_ncvmp.R")

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
                           data = riesby)
summary(riesby_ncvmp)

riesby_ncvmp_boot <- bootstrap_mels_ncvmp(riesby_ncvmp, B = 1000)
summary(riesby_ncvmp_boot)

###########################################
### Example 1. Health Behavior Data Set ###
###########################################

indat <- read.table("Dataset_HealthBehavAcadPerfAffect.dat", header = FALSE, 
                    col.names=c("id", "day", "sex", "age", "sem", "sq", "physact", "pa", "na", "lga", 
                                "exam", "hsg", "bdi", "day_c"), na.strings="-99")

indat_ncvmp <- mels_ncvmp(y = "pa", 
                          beta_formula = ~ day_c, 
                          alpha_formula = ~ sex, 
                          tau_formula = ~ day_c, 
                          id = "id", 
                          data = indat)
summary(indat_ncvmp)

indat_ncvmp_boot <- bootstrap_mels_ncvmp(indat_ncvmp, B = 1000)
summary(indat_ncvmp_boot)

posmod <- read.table("posmood_example.dat", header = FALSE, 
                     col.names = c("id", "posmood", "alone", "genderf"))
posmod_ncvmp <- mels_ncvmp(y = "posmood", 
                           beta_formula = ~ alone + genderf, 
                           alpha_formula = ~ genderf, 
                           tau_formula = ~ alone + genderf, 
                           id = "id", 
                           data = posmod)
summary(posmod_ncvmp)

