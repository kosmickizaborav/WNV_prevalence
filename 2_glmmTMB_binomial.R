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

wnv_dt <- readRDS(file.path(cluster_dir, "1_WNV_prevalence_data_model.rds"))

wnv_dt[, avilist_genus := sub(" .*", "", avilist_name)]
wnv_dt <- wnv_dt[, negative := total_tested - positive]

setwd(model_dir)


add_continents <- function(
    steps, crs = sf::st_crs(4326), align_start = T,
    coord_cols = NULL, scale = "medium"){

  # overlay the map with start or the end point of the step (didn't check both
  # to preserve the country information in case we need it later)
  if(is.null(coord_cols)){
    coord_cols <- if(align_start) c("x1_", "y1_") else c("x2_", "y2_")
  }

  steps <- steps[
    , geometry := sf::st_as_sf(
      .SD, coords = coord_cols, crs = crs), .SDcols = coord_cols]

  # load or download the continent data with the right resolution
  conti <- tryCatch(
    rnaturalearth::ne_load(
      type = "geography_regions_polys",  scale = scale, category = "physical"),
    error = function(e) {
      rnaturalearth::ne_download(
        type = "geography_regions_polys",  scale = scale, category = "physical")
    }
  )

  # check the nearest geometry
  conti_id <- sf::st_nearest_feature(steps$geometry, conti)

  # convert to datatable for merging
  conti_dt <- as.data.table(conti)[
    , .(continent = REGION, subregion = SUBREGION)]
  conti_dt <- conti_dt[conti_id]

  # add columns to the original step data
  steps <- cbind(steps[, geometry := NULL], conti_dt)

  return(steps)

}


wnv_dt <- add_continents(
  wnv_dt, coord_cols = c("long", "lat"), scale = "medium")



# Frequentist -------------------------------------------------------------

fixed_eff <- c(
  "1", "tarsus_log", "freshwater", "primary_lifestyle", 
   "trophic_niche", "migration", "sociality", "clutch_max_log", "nest_placement", 
  "longevity_log", "abundance_log", "log_human_population_density"
  )

base_f <- "cbind(positive, negative) ~ %s (1 | avilist_order/avilist_family/avilist_name) + (1 | method_cat) + (1 | continent/country) + (1 | sampling_year)" 

for(i in seq_along(fixed_eff)){
  
  all_fe <- paste(fixed_eff[1:i], collapse = " + ")
  
  formulaf <- sprintf(base_f, paste(all_fe, "+"))
  
  mod <- glmmTMB(
    as.formula(formulaf),
    data = wnv_dt,
    family = binomial(),
    ziformula = ~1
  )
  
  mname <- paste0("2_", i, paste(fixed_eff[1:i], collapse = "_"), ".rds")
  
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
# 7: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 8: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 9: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 10: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 11: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 12: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 13: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 14: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 15: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 16: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 17: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 18: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 19: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 20: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 21: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 22: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 23: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 24: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 25: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 26: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 27: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 28: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 29: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 30: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 31: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 32: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 33: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 34: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 35: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 36: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 37: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 38: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 39: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 40: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 41: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 42: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
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
# 48: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 49: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation
# 50: In (function (start, objective, gradient = NULL, hessian = NULL,  ... :
#   NA/NaN function evaluation



# Full formula ------------------------------------------------------------

full_f <- "positive | trials(total_tested) ~ 1 + tarsus_log + longevity_log + clutch_max_log + log_human_population_density + abundance_log + freshwater + migration + sociality + primary_lifestyle + nest_placement + trophic_niche + (1 | gr(avilist_name,cov=A)) + (1 | method_cat) + (1 | gr(location_id, cov = Sigma_spatialAF)) + (1 | sampling_year)"

m1 <- glmmTMB(
  cbind(positive, total_tested-positive) ~ tarsus_log + longevity_log + clutch_max_log + log_human_population_density + abundance_log + freshwater + migration + sociality + primary_lifestyle + nest_placement + trophic_niche + (1 | avilist_name) + (1 | avilist_order) + (1 | method_cat) + (1 | country:sampling_year),
  data = wnv_dt,
  family = binomial(),
  ziformula = ~1
)


# check models ------------------------------------------------------------

files <- list.files(model_dir, full.names = T)

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
AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12)
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
