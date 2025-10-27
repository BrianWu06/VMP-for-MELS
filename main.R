library(lme4)
library(foreach)
library(doParallel)
source("mels_ncvmp.R")
source("mels_ncvmp_old.R")

##################################
### Example 1. Riesby Data Set ###
##################################

riesby <- read.table("RIESBY.DAT.txt", na.strings = ".", 
                     col.names = c("id", "hamd", "intcpt", "week", "endog", "endweek"))

riesby_ncvmp <- mels_ncvmp(y = "hamd", 
                           beta_formula = ~ week + endog + endweek, 
                           alpha_formula = ~ endog, 
                           tau_formula = ~ week + endog, 
                           id = "id",
                           data = riesby)
summary(riesby_ncvmp)

riesby_ncvmp_boot <- bootstrap_mels_ncvmp(riesby_ncvmp, B = 1000, cores = 10)
summary(riesby_ncvmp_boot)

###########################################
### Example 2. Health Behavior Data Set ###
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

indat_ncvmp_boot <- bootstrap_mels_ncvmp(indat_ncvmp, B = 1000, cores = 10)
summary(indat_ncvmp_boot)


###########################################
### Example 3. Positive Mood Data Set #####
###########################################

posmood <- read.table("moods_example_augmented.dat", header = FALSE, 
                       col.names = c("id", "posmood", "negmood", "t1", "t2", "t3", "t4", "w1", 
                                     "w2", "w3", "w4", "w5", "w6", "other_bs", "other_ws", 
                                      "genderf", "age15", "tirbor", "frustr"))

posmood_ncvmp <- mels_ncvmp(y = "posmood", 
                            beta_formula = ~ other_bs + other_ws + genderf + t1 + t2 +
                              t3 + t4 + w1 + w2 + w3 + w4 + w5 + w6 + tirbor + frustr, 
                            alpha_formula = ~ other_bs + genderf + age15, 
                            tau_formula = ~ other_bs + other_ws + genderf + age15 + tirbor + frustr, 
                            id = "id", 
                            data = posmood)
summary(posmood_ncvmp)


posmood_ncvmp_boot <- bootstrap_mels_ncvmp(posmood_ncvmp, B = 1000, cores = 10)
summary(posmood_ncvmp_boot)
