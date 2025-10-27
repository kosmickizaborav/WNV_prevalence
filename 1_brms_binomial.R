library(data.table)
library(ape)
library(brms)

laptop <- F
random <- F

main_dir <- "/home/nbogdanovic/WNV_prevalence"
cluster_dir <- file.path(main_dir, "Data_for_cluster")
model_dir <- file.path(main_dir, "Models")


if(laptop == T){
  cluster_dir <- here::here("Data_for_cluster")
  model_dir <- here::here("Models")
}

dir.create(model_dir, showWarnings = F)

# Seed stuff
set.seed(202510)
BAYES_SEED <- 202510

# 0: Prepare data PC---------------------------------------------------------

# # laptop
# data_dir <- here::here("Data")
# # data for cluster dir
# cluster_dir <- here::here("Data_for_cluster")
# dir.create(cluster_dir, showWarnings = F)
# 
# wnv_dt <- fread(file.path(data_dir, "00_WNV_prevalence_data.csv"))[
#   , name_type := NULL]
# 
# wnv_dt <- wnv_dt[
#   , ":=" (
#     avilist_name = gsub(" ", "_", avilist_name),
#     method = gsub(" ", "_", method)
#   )]
# 
# # selected traits
# traits_select <- c(
#   "birdlife_name", "mass_log", "tarsus_log", "hwi_log", "longevity_log",
#   "clutch_max_log", "log_human_population_density", "abundance_log",
#   "habitat", "habitat_TP2019", "freshwater", "migration", "sociality",
#   "primary_lifestyle", "nest_placement", "trophic_niche", "altitude_cat"
# )
# 
# traits_dt <- fread(file.path(data_dir, "00_traits_filtered.csv"))
# 
# traits_dt <- traits_dt[
#   , ':=' (
#     mass_log = scale(log(mass)),
#     tarsus_log = scale(log(tarsus_length)),
#     hwi_log = scale(log(hand_wing_index)),
#     longevity_log = scale(log(maximum_longevity)),
#     clutch_max_log = scale(log(clutch_max)),
#     abundance_log = scale(log10(abundance_estimate + 1)),
#     log_human_population_density = scale(log_human_population_density),
#     migration = factor(
#       fifelse(migration == "sedentary", "sedentary", "migratory")),
#     nest_placement = factor(nest_placement),
#     habitat = factor(fcase(
#       habitat == "Human Modified", "humanModified",
#       habitat %in% c("Forest", "Woodland", "Shrubland", "Grassland"), "greenLandscape",
#       habitat %in% c("Wetland", "Riverine", "Coastal"), "waterLandscape",
#       habitat %in% c("Desert", "Rock"), "rockyLandscape",
#       habitat == "Marine", "marineLandscape"
#     )),
#     sociality = factor(fcase(
#       colonial == 1 | social == 1, "social",
#       default = "nonSocial"
#     )),
#     primary_lifestyle = factor(fcase(
#       primary_lifestyle %in% c("Aerial", "Insessorial"), "air",
#       primary_lifestyle == "Aquatic", "water",
#       default = "other"
#     )),
#     trophic_niche = factor(fcase(
#       trophic_niche %in% c("Scavenger", "Vertivore"), "carnivore",
#       trophic_niche %in% c("Aquatic predator",  "Herbivore aquatic"), "aquaticFeeder",
#       trophic_niche %in% c("Frugivore", "Granivore", "Nectarivore"), "cropFeeder",
#       default = "terrestrialFeeder")),
#     altitude_cat = factor(fcase(
#       minimum_altitude < 500, "lowland",
#       minimum_altitude >= 500 & minimum_altitude < 1500, "midland",
#       minimum_altitude >= 1500, "highland"
#     ))
#   )][, ..traits_select]
# 
# wnv_dt <- merge(wnv_dt, traits_dt, by = "birdlife_name")
# 
# 
# # keep only the data that is complete
# wnv_dt[
#   , data_complete := apply(.SD, 1, function(row) all(!is.na(row) & row != ""))]
# wnv_dt <- wnv_dt[data_complete == T]
# 
# # and minimum 20 tested individuals per species
# wnv_dt[, ind_per_sp := sum(total_tested), by = avilist_name]
# wnv_dt <- wnv_dt[ind_per_sp >= 20]
# 
# 
# big_tests <- c("serological_ELISA", "molecular_detection", "serological_VNT")
# wnv_dt[, method_cat := fifelse(!method %in% big_tests, "other", method)]
# 
# # save data
# fwrite(wnv_dt, file.path(cluster_dir, "1_WNV_prevalence_data_model.csv"))

# #PHYLOGENETIC AUTOCORRELATION
# bird_tree <- readRDS(file.path(data_dir, "00_bird_tree_match_avilist.rds"))
# bird_tree <- drop.tip(
#   bird_tree,
#   bird_tree$tip.label[!bird_tree$tip.label %in% unique(wnv_dt$avilist_name)]
#   )
# # this way we generate covariance matrix that will be inputted in the model as:
# # gr(scientific_name_bt, cov = A) - we make sure that species are correlated
# # as specified by the covariance matrix A
# A <- vcv.phylo(bird_tree)
# saveRDS(A, file.path(cluster_dir, "1_bird_tree_A.rds"))


# TRIED THIS, BUT THE DISTANCES HAVE BIMODAL DISTRIBUION,
# it seems that this spaial autocorrelation matric is not a good idea here
# skipped it, for now
#
# study_locs <- unique(wnv_dt[, .(country, long, lat)])
# setorder(study_locs, country, long, lat)
# study_locs[
#   , location_id := paste0(gsub(" ", "_", country), "_", 1:.N), by = country]
# study_locs <- sf::st_as_sf(study_locs, coords = c("long", "lat"), crs = 4326)
# spatiald <- as.numeric(sf::st_distance(study_locs))/1000 # in km
# dist_matrix_mAF <- distm(
#   x = unique_coordsAF[, c("decimalLongitude", "decimalLatitude")],  # lon, lat
#   fun = distHaversine) # in M
# # in KM
# dist_matrix_kmAF <- dist_matrix_mAF / 1000
# # Rename matrix with location_id
# rownames(dist_matrix_kmAF) <- unique_coordsAF$location_id
# colnames(dist_matrix_kmAF) <- unique_coordsAF$location_id
# range_paramAF <- 500  # decay scale in km
# Sigma_spatialAF <- exp(-dist_matrix_kmAF / range_paramAF)


# add_continents <- function(
#     steps, crs = sf::st_crs(4326), align_start = T,
#     coord_cols = NULL, scale = "medium"){
#
#   # overlay the map with start or the end point of the step (didn't check both
#   # to preserve the country information in case we need it later)
#   if(is.null(coord_cols)){
#     coord_cols <- if(align_start) c("x1_", "y1_") else c("x2_", "y2_")
#   }
#
#   steps <- steps[
#     , geometry := sf::st_as_sf(
#       .SD, coords = coord_cols, crs = crs), .SDcols = coord_cols]
#
#   # load or download the continent data with the right resolution
#   conti <- tryCatch(
#     rnaturalearth::ne_load(
#       type = "geography_regions_polys",  scale = scale, category = "physical"),
#     error = function(e) {
#       rnaturalearth::ne_download(
#         type = "geography_regions_polys",  scale = scale, category = "physical")
#     }
#   )
#
#   # check the nearest geometry
#   conti_id <- sf::st_nearest_feature(steps$geometry, conti)
#
#   # convert to datatable for merging
#   conti_dt <- as.data.table(conti)[
#     , .(continent = REGION, subregion = SUBREGION)]
#   conti_dt <- conti_dt[conti_id]
#
#   # add columns to the original step data
#   steps <- cbind(steps[, geometry := NULL], conti_dt)
#
#   return(steps)
#
# }
#
#
# dist_locs <- add_continents(
#   wnv_dt, coord_cols = c("long", "lat"), scale = "medium")



# 0: Load data cluster -------------------------------------------------------

A <- readRDS(file.path(cluster_dir, "1_bird_tree_A.rds"))

wnv_dt <- fread(file.path(cluster_dir, "1_WNV_prevalence_data_model.csv"))

setwd(model_dir)

# 1: Intercept + random effects-------------------------------------------------

if(random == T){
  
  # get the log warning files
  logfile <- file.path(model_dir, "all_model_warnings.log")
  logcon <- file(logfile, open = "wt")
  
  # List the random effects you want to test
  random_effects <- c(
    "(1 | gr(avilist_name,cov=A))",
    "(1 | method_cat)",
    "(1 | country:sampling_year)"
  )

  for(re in random_effects) {

    var_id <- gsub("[^[:alnum:]]", "_", gsub("^\\(1 \\| |\\)$", "", re))

    mname <- paste0("1_binomial_r_", var_id, ".rds")

    # Dynamically build the formula string
    formula_str <- paste0("positive | trials(total_tested) ~ 1 + ", re)

    if(re != "(1 | gr(avilist_name,cov=A))"){
      formula_str <- paste0(
        "positive | trials(total_tested) ~ 1 + (1 | gr(avilist_name,cov=A)) + ",
        re)
    }

    # convert to brms formula object
    f <- bf(formula(formula_str), family = binomial())

    # my original priors
    # p <- c(
    #   prior(normal(0, 1.5), class = Intercept),
    #   prior(exponential(1), class = sd)
    # )

    # uninformative priors, picked from the comparative_analyais_brm_clean
    p = c(
      prior(student_t(3, 0, 10), "Intercept"),
      prior(student_t(3, 0, 10), "sd")
    )

    # p <- get_prior(f, data = wnv_dt, data2 = list(A = A))

    # Fit the model only if not already saved
    if (!file.exists(mname)) {
      withCallingHandlers({
        m <- brm(
          data = wnv_dt,
          data2 = list(A = A),
          formula = f,
          prior = p,
          iter = 6000,
          warmup = 2000,
          cores = 4,
          chains = 4,
          seed = BAYES_SEED,
          sample_prior = TRUE,
          file = mname
        )
      }, warning = function(w) {
        writeLines(
          paste0(Sys.time(), " | MODEL: ", mname, " \n ", conditionMessage(w)),
          logcon
        )
        invokeRestart("muffleWarning")
      })
    }

    writeLines(c(" ", "MODEL:", mname, "DONE!", " "))

    # TRY
    mname <- paste0("1_zero-negbinomial_r_", var_id, ".rds")

    # Dynamically build the formula string
    formula_str <- paste0("positive ~ 1 + total_tested + ", re)

    if(re != "(1 | gr(avilist_name,cov=A))"){
      formula_str <- paste0(
        "positive ~ 1 + total_tested + (1 | gr(avilist_name,cov=A)) + ", re)
    }

    # convert to brms formula object
    f <- bf(formula(formula_str), family = zero_inflated_negbinomial())

    # I ran the script with this called prior instead of p, but lickly the 
    # shape and zi parameters are the same in default priors so shouldn't 
    # matter
    p <- c(
      prior(student_t(3, 0, 10), class = "Intercept"),
      prior(student_t(3, 0, 10), class = "sd"),
      prior(gamma(0.01, 0.01), class = "shape"),      # for negative binomial shape
      prior(beta(1,1), class="zi")                    # only for zero-inflated
    )

    # Fit the model only if not already saved
    if (!file.exists(mname)) {
      withCallingHandlers({
        m <- brm(
          data = wnv_dt,
          data2 = list(A = A),
          formula = f,
          prior = p,
          iter = 6000,
          warmup = 2000,
          cores = 4,
          chains = 4,
          seed = BAYES_SEED,
          sample_prior = TRUE,
          file = mname
        )
      }, warning = function(w) {
        writeLines(
          paste0(Sys.time(), " | MODEL: ", mname, " \n ", conditionMessage(w)),
          logcon
        )
        invokeRestart("muffleWarning")
      })
    }

    writeLines(c(" ", "MODEL:", mname, "DONE!", " "))

  }


  close(logcon) # Close log file connection

  
}


# 2: Intercept + fixed effects-------------------------------------------------

setwd(model_dir)

# get the log warning files
logfile <- file.path(model_dir, "all_model_fixed_warnings.log")
logcon <- file(logfile, open = "wt")


# selected traits
traits_select <- c(
  "birdlife_name", "mass_log", "tarsus_log", "longevity_log",
  "clutch_max_log", "log_human_population_density", "abundance_log",
  "habitat", "freshwater", "migration", "sociality",
  "primary_lifestyle", "nest_placement", "trophic_niche"
)

body_size <- c("mass_log", "tarsus_log", "longevity_log")
habitat_vars <- c(
  "habitat", "freshwater", "primary_lifestyle", "trophic_niche")

other_vars <- setdiff(
  traits_select, c(body_size, habitat_vars, "birdlife_name"))

# Generate all possible subsets of other_vars (including empty subset)
other_combos <- unlist(lapply(0:length(other_vars), function(k) {
  combn(other_vars, k, simplify = FALSE)
}), recursive = FALSE)

# Build formulas: each must have one body_size, one habitat_var, plus any other_vars
formulas <- c()
i <- 1
for (bs in body_size) {
  for (hv in habitat_vars) {
    for (ov in other_combos) {
      model_vars <- c(bs, hv, ov)
      formulas[i] <- paste(model_vars, collapse = " + ")
      i <- i + 1
    }
  }
}

length(formulas)
formulas_dt <- data.table(formula = formulas)[, fname := paste0("f", 1:.N)]

fwrite(formulas_dt, "1_model_formulas_fixed_effects.csv")


re <- "(1 | gr(avilist_name,cov=A)) + (1 | method_cat) + (1 | country:sampling_year)"

# TRY
for(i in 1:nrow(formulas_dt)){
  
  fn <- formulas_dt[i, fname]
  fp <- formulas_dt[i, formula]
  
  mname <- paste0("1_zero-negbinomial_f_", fn, ".rds")
  
  # Dynamically build the formula string
  formula_str <- paste0("positive ~ 1 + total_tested + ", fp, " + ", re)
  
  # convert to brms formula object
  f <- bf(formula(formula_str), family = zero_inflated_negbinomial())
  
  p <- c(
    prior(student_t(3, 0, 10), class = "Intercept"),
    prior(student_t(3, 0, 10), class = "sd"),
    # for negative binomial shape
    prior(gamma(0.01, 0.01), class = "shape"),
    prior(beta(1,1), class="zi")
    # originally prior here was not specified so I just left it like it is
    # prior(normal(0, 2), class = "b")
  )
  
  # Fit the model only if not already saved
  if (!file.exists(mname)) {
    withCallingHandlers({
      m <- brm(
        data = wnv_dt,
        data2 = list(A = A),
        formula = f,
        prior = p,
        iter = 6000,
        warmup = 2000,
        cores = 8,
        chains = 4,
        seed = BAYES_SEED,
        sample_prior = TRUE,
        file = mname
      )
    }, warning = function(w) {
      writeLines(
        paste0(Sys.time(), " | MODEL: ", mname, " \n ", conditionMessage(w)),
        logcon
      )
      invokeRestart("muffleWarning")
    })
  }
  
  writeLines(c(" ", "MODEL:", mname, "DONE!", " "))
  
}


close(logcon) # Close log file connection

  
  




# check models ------------------------------------------------------------

# models_dir <- here::here("Models_first_trial")
# 
# 
# mnames <- list.files(models_dir, pattern = ".rds", full.names = T)
# all_models <- lapply(mnames, readRDS)
# 
# 
# summarym <- lapply(all_models, summary)
# loom <- lapply(all_models, loo)
# 
# loo_compare(loom)
# 
# summary(m)
 
# # variables ---------------------------------------------------------------
# 
# categorical_vars <- c("migration", "trophic_niche",  "primary_lifestyle",
#                       "nest_placement", "foraging", "habitat")
# 
# 
# # 2: Variables steps-----------------------------------------------------------
# 
# for(var in categorical_vars){
#   
#   bird_db_var <- bird_db |> 
#     select(all_of(c("bt_name", "total", "positive", "philo", var))) |> 
#     rename_with(~"vari", all_of(var))
#   
#   mname <- str_c("01_mph_1_", var)
#   
#   if(!file.exists(str_c(mname, ".rds"))) {
#     
#     f <- bf(
#       positive | trials(total) ~ a + b, 
#       a ~ 1 + (1 | bt_name) + (1 | gr(philo, cov = A)),
#       b ~ 0 + vari, 
#       nl = TRUE
#     )
#     
#     p <- c(
#       prior(normal(0, 0.5), nlpar = b),
#       prior(normal(0, 1.5), class = b, coef = Intercept, nlpar = a),
#       prior(exponential(1), class = sd, nlpar = a)
#     )
#     
#     m <- brm(
#       data = bird_db_var,
#       data2 = list(A = A),
#       family = binomial,
#       formula = f, 
#       prior = p, 
#       iter = 5000, # iterations including the burn-in
#       warmup = 2000, # warm-up, the burn-in period,
#       cores = 4,
#       chains = 4,
#       seed = BAYES_SEED,
#       sample_prior = T,
#       # control = list(adapt_delta = .99),
#       # backend = "cmdstanr",
#       file = mname
#     )
#     
#     warnings()
#     
#     rm(m)
#     
#     writeLines(c(" ", "MODEL:", mname, "DONE!", " ")) 
#     
#     
#     
#   }
#   
# }
# 
# # Warning messages:
# #   1: There were 12000 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# # https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# # 2: Examine the pairs() plot to diagnose sampling problems
# 
# 
# # Partitioning categorical variance ------------------------------------------
# 
# # ADDED TO THE SCRIPT: 10/04/2024
# 
# #get_prior(f, data = bird_db, data2 = list(A = A), family = binomial)
# 
# mname <- "01_mph_1_catvar_rand"
# 
# f <- bf(
#   positive | trials(total) ~ 1 + (1|habitat) + (1|migration) + (1|trophic_niche) + (1|primary_lifestyle) + (1|nest_placement) + (1|foraging) + (1 | bt_name) + (1 | gr(philo, cov = A))
# )
# 
# p <- c(
#   prior(normal(0, 1.5), class = Intercept),
#   prior(exponential(1), class = sd)
# )
# 
# m <- brm(
#   data = bird_db,
#   data2 = list(A = A),
#   family = binomial,
#   formula = f, 
#   prior = p, 
#   iter = 5000, # iterations including the burn-in
#   warmup = 2000, # warm-up, the burn-in period,
#   cores = 4,
#   chains = 4,
#   seed = BAYES_SEED,
#   sample_prior = T,
#   # control = list(adapt_delta = .99),
#   # backend = "cmdstanr",
#   file = mname
# )
# 
# writeLines(c(" ", "MODEL:", mname, "DONE!", " ")) 
# 
# warnings()
# 
