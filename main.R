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

# riesby_clean <- na.omit(riesby[, c("hamd", "week", "endog", "id", "endweek")])
# riesby_y <- riesby_clean$hamd
# riesby_x <- model.matrix(~ week + endog + endweek, data = riesby_clean)
# riesby_z <- model.matrix(~ week + endog, data = riesby_clean)
# riesby_id <- as.numeric(factor(riesby_clean$id))
# riesby_fastregls <- fastregls(riesby_y, riesby_x, riesby_z, riesby_id, getSigmaSE = FALSE)
# riesby_fastregls$beta
# riesby_fastregls$tau
# riesby_fastregls$sigma_sq_nu
# riesby_fastregls$sigma_sq_omega

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
                          r = FALSE, 
                          data = indat)
summary(indat_ncvmp)

indat_ncvmp_r <- mels_ncvmp(y = "pa", 
                            beta_formula = ~ day_c, 
                            alpha_formula = ~ sex, 
                            tau_formula = ~ day_c, 
                            id = "id", 
                            r = TRUE, 
                            data = indat)
summary(indat_ncvmp_r)

indat_ncvmp_boot <- bootstrap_mels_ncvmp(indat_ncvmp, B = 1000)
summary(indat_ncvmp_boot)


###########################################
### Example 3. Positive Mood Data Set #####
###########################################

posmood <- read.table("moods_example_augmented.dat", header = FALSE, 
                       col.names = c("id", "posmood", "negmood", "t1", "t2", "t3", "t4", "w1", 
                                     "w2", "w3", "w4", "w5", "w6", "other_bs", "other_ws", 
                                      "genderf", "age15", "tirbor", "frustr"))
posmood_ncvmp <- mels_ncvmp(y = "posmood", 
                            beta_formula = ~ other_bs + other_ws + genderf + age15 + tirbor + frustr, 
                            alpha_formula = ~ genderf + age15, 
                            tau_formula = ~ other_bs + other_ws + genderf + age15, 
                            r = FALSE, 
                            id = "id", 
                            data = posmood)
summary(posmood_ncvmp)

posmood_ncvmp_r <- mels_ncvmp(y = "posmood", 
                              beta_formula = ~ other_bs + other_ws + genderf + age15 + tirbor + frustr, 
                              alpha_formula = ~ genderf + age15, 
                              tau_formula = ~ other_bs + other_ws + genderf + age15, 
                              r = TRUE, 
                              id = "id", 
                              data = posmood)
summary(posmood_ncvmp_r)
