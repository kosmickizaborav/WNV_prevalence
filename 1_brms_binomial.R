library(data.table)
library(ape)
library(brms)

data_dir <- here::here("Data")
model_dir <- file.path(data_dir, "Models")
dir.create(model_dir, showWarnings = F)

# Seed stuff
set.seed(202510)
BAYES_SEED <- 202510

# 0: Prepare data ---------------------------------------------------------

# # laptop
# data_dir <- here::here("Data")
# # data for cluster dir
# cluster_dir <- here::here("Data_for_cluster")
# dir.create(cluster_dir, showWarnings = F)
# 
# wnv_dt <- fread(file.path(data_dir, "00_WNV_prevalence_data.csv"))
# wnv_dt[, ind_per_sp := sum(total_tested), by = avilist_name]
# 
# wnv_dt <- wnv_dt[ind_per_sp >= 20][
#   , avilist_name := gsub(" ", "_", avilist_name)]
# 
# # selected traits
# traits_select <- c(
#   "birdlife_name", "mass_log", "tarsus_log", "hwi_log", "longevity_log",
#   "clutch_max_log", "log_human_population_density",
#   "habitat", "habitat_TP2019", "migration", "freshwater",
#   "primary_lifestyle", "nest_placement", "trophic_niche"
# )
# 
# traits_dt <- fread(file.path(data_dir, "00_traits_filtered.csv"))
# 
# traits_dt <- traits_dt[
#   , ':=' (
#     mass_log = log(mass),
#     tarsus_log = log(tarsus_length),
#     longevity_log = log(maximum_longevity),
#     clutch_max_log = log(clutch_max),
#     hwi_log = log(hand_wing_index), 
#     abundance_log10 = log10(abundance_estimate), 
#     migration = factor(
#       fifelse(migration == "sedentary", "sedentary", "migratory")),
#     nest_placement = factor(nest_placement),
#     habitat = factor(fcase(
#       habitat == "Human Modified", "humanModified", 
#       habitat %in% c("Woodland", "Shrubland", "Grassland"), "greenLandscape",
#       habitat %in% c("Wetland", "Riverine", "Coastal"), "waterLandscape",
#       habitat %in% c("Desert", "Rocky"), "rockyLandscape",
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
#   )][, ..traits_select][birdlife_name %in% unique(wnv_dt$birdlife_name)]
# 
# traits_dt[
#   , data_complete := apply(.SD, 1, function(row) all(!is.na(row) & row != ""))]
# 
# traits_dt <- traits_dt[data_complete == T]
# 
# wnv_dt <- wnv_dt[birdlife_name %in% traits_dt$birdlife_name]
# 
# wnv_dt <- merge(wnv_dt, traits_dt, by = "birdlife_name")
# fwrite(wnv_dt, file.path(cluster_dir, "1_WNV_prevalence_data_model.csv"))

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


# PHYLOGENETIC AUTOCORRELATION
bird_tree <- readRDS(file.path(data_dir, "00_bird_tree_match_avilist.rds"))
bird_tree <- ape::drop.tip(
  bird_tree, 
  bird_tree$tip.label[!bird_tree$tip.label %in% unique(wnv_dt$avilist_name)]
  )
# this way we generate covariance matrix that will be inputted in the model as:
# gr(scientific_name_bt, cov = A) - we make sure that species are correlated
# as specified by the covariance matrix A
A <- vcv.phylo(bird_tree)
saveRDS(A, file.path(cluster_dir, "1_bird_tree_A.rds"))


# 1: Intercept + random effects-------------------------------------------------

# List the random effects you want to test
random_effects <- c(
  "(1 | gr(avilist_name,cov=A))",
  "(1 | country)",
  "(1 | sampling_year)",
  "(1 | method)", 
  "(1 | country:sampling_year)"
)

for(re in random_effects) {
  
  mname <- paste0(
    "1_r_", gsub("[^[:alnum:]]", "_", gsub("^\\(1 \\| |\\)$", "", re)), ".rds")
  
  # Dynamically build the formula string
  formula_str <- paste0("positive | trials(total_tested) ~ 1 + ", re)
  
  if(re != "(1 | gr(avilist_name,cov=A))"){
    formula_str <- paste0(
      "positive | trials(total_tested) ~ 1 + (1 | gr(avilist_name,cov=A)) + ",
      re)
  }
  
  # convert to brms formula object
  f <- bf(as.formula(formula_str))
  
  # Fit the model only if not already saved
  if (!file.exists(mname)) {
    
    m <- brm(
      data = wnv_dt,
      data2 = list(A = A),
      family = binomial,
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
  }
  
  writeLines(c(" ", "MODEL:", mname, "DONE!", " "))

  warnings()

  
}


 
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
