library(data.table)
library(DHARMa)
library(glmmTMB)

prepare_data <- F

# here package should read it automatically from R project so no need to specify
# the base directory
#main_dir <- "/home/nbogdanovic/WNV_prevalence"
# cluster_dir <- file.path(main_dir, "Data_for_cluster")
# model_dir <- file.path(main_dir, "Models")

cluster_dir <- here::here("Data_for_cluster")
model_dir <- here::here("Models_freq")
dir.create(model_dir, showWarnings = F)

# 0: Load data cluster -------------------------------------------------------

A <- readRDS(file.path(cluster_dir, "1_bird_tree_A.rds"))

wnv_dt <- readRDS(file.path(cluster_dir, "1_WNV_prevalence_data_model.rds"))

setwd(model_dir)



# Frequentist -------------------------------------------------------------

fixed_eff <- c(
  "", "tarsus_log", "freshwater", "migration", "sociality", 
  "nest_placement",  "clutch_max_log",  "longevity_log", "trophic_niche",  
  "abundance_log", "log_human_population_density", "altitude_cat",  "primary_lifestyle"
  )


wnv_dt <- wnv_dt[, negative := total_tested - positive]

base_f <- "cbind(positive, negative) ~ 1 + %s(1 | avilist_family) + (1 | avilist_order) + (1 | method_cat) + (1 | country:sampling_year)" 

for(i in seq_along(fixed_eff)){
  
  fe <- fixed_eff[i]
  if(i > 1){
    all_fe <- paste(fixed_eff[2:i], collapse = " + ")
  }
  
  formulaf <- sprintf(base_f, ifelse(fe == "", "", paste(all_fe, "+")))
  
  mod <- glmmTMB(
    as.formula(formulaf),
    data = wnv_dt,
    family = binomial(),
    ziformula = ~1
  )
  
  mname <- paste0(
    "2_2_", ifelse(i==1, "intercept", paste(fixed_eff[2:i], collapse = "_")), ".rds")
  
  saveRDS(mod, file.path(model_dir, mname))
  
}



# > warnings()
# Warning messages:
# 1: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 2: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 3: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 4: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 5: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 6: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 7: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 8: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 9: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 10: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 11: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 12: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 13: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 14: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 15: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 16: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 17: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 18: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 19: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 20: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 21: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 22: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 23: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 24: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 25: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 26: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 27: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 28: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 29: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 30: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 31: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 32: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 33: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 34: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 35: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 36: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 37: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 38: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 39: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 40: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 41: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 42: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 43: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 44: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 45: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 46: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 47: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 48: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 49: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 50: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation


# check models ------------------------------------------------------------

files <- list.files(model_dir, "^2_2_tarsus_", full.names = T)

all_models <- lapply(files, readRDS)

m1 <- all_models[[1]]
m2 <- all_models[[2]]
m3 <- all_models[[3]]
m4 <- all_models[[4]]
m5 <- all_models[[5]]
m6 <- all_models[[6]]
m7 <- all_models[[7]]
m8 <- all_models[[8]]
m9 <- all_models[[9]]
m10 <- all_models[[10]]
m11 <- all_models[[11]]
m12 <- all_models[[12]]

anova(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, test = "Chisq")
AIC(unlist(all_models))
# Compare models using AIC
AIC(m1, m2, m3, m4, m5, m6)

# Compare models using BIC
BIC(m1, m2, m3)


anova(m1, m2, m3, m4, test = "Chisq") 


simres <- simulateResiduals(binbinfit)
plot(simres)  # Residual plots, tests for overdispersion, outliers, zero-inflation


simresz <- simulateResiduals(fit_zi)
plot(simresz) 

# For probabilities:
pred_probs <- predict(fit_zi, type = "response")  # This gives you predicted probabilities

# For predicted counts (optional, for binomial-type models):
pred_counts <- pred_probs * wnv_dt$total_tested

plot(
  wnv_dt$positive,   # x: observed counts
  pred_counts,          # y: predicted counts
  xlab = "Observed Counts",
  ylab = "Predicted Counts",
  main = "Predicted vs Observed Counts"
)
abline(0, 1, col = "red")
