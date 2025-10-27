library(data.table)
library(ape)
library(brms)

prepare_data <- F

# here package should read it automatically from R project so no need to specify
# the base directory
#main_dir <- "/home/nbogdanovic/WNV_prevalence"
# cluster_dir <- file.path(main_dir, "Data_for_cluster")
# model_dir <- file.path(main_dir, "Models")

cluster_dir <- here::here("Data_for_cluster")
model_dir <- here::here("Models")

dir.create(model_dir, showWarnings = F)

# Seed stuff
set.seed(202510)
BAYES_SEED <- 202510

# 0: Prepare data PC---------------------------------------------------------

if(prepare_data == T){
  
  # # laptop
  data_dir <- here::here("Data")
  # data for cluster dir
  cluster_dir <- here::here("Data_for_cluster")
  dir.create(cluster_dir, showWarnings = F)
  
  wnv_dt <- fread(file.path(data_dir, "00_WNV_prevalence_data.csv"))[
    , name_type := NULL]
  
  wnv_dt <- wnv_dt[
    , ":=" (
      avilist_name = gsub(" ", "_", avilist_name),
      method = gsub(" ", "_", method)
    )]
  
  # selected traits
  # removed "habitat_TP2019",
  traits_select <- c(
    "birdlife_name", "mass_log", "tarsus_log", "hwi_log", "longevity_log",
    "clutch_max_log", "log_human_population_density", "abundance_log",
    "habitat", "freshwater", "migration", "sociality",
    "primary_lifestyle", "nest_placement", "trophic_niche", "altitude_cat"
  )
  
  traits_dt <- fread(file.path(data_dir, "00_traits_filtered.csv"))
  
  traits_dt <- traits_dt[
    , ':=' (
      mass_log = scale(log(mass)),
      tarsus_log = scale(log(tarsus_length)),
      hwi_log = scale(log(hand_wing_index)),
      longevity_log = scale(log(maximum_longevity)),
      clutch_max_log = scale(log(clutch_max)),
      abundance_log = scale(log(abundance_estimate + 1)),
      log_human_population_density = scale(log_human_population_density),
      migration = factor(
        fifelse(migration == "sedentary", "sedentary", "migratory")),
      nest_placement = factor(nest_placement),
      habitat = fcase(
        habitat == "Human Modified", "humanModified",
        habitat %in% c("Forest", "Woodland", "Shrubland", "Grassland"), "greenLandscape",
        habitat %in% c("Wetland", "Riverine", "Coastal"), "waterLandscape",
        habitat %in% c("Desert", "Rock"), "rockyLandscape",
        habitat == "Marine", "marineLandscape"
      ),
      freshwater = fifelse(freshwater == 1, "freshwater", "other"),
      sociality = fcase(
        colonial == 1 | social == 1, "social",
        default = "nonSocial"),
      primary_lifestyle = fcase(
        primary_lifestyle %in% c("Aerial", "Insessorial"), "air",
        primary_lifestyle == "Aquatic", "water",
        default = "other"
      ),
      trophic_niche = fcase(
        trophic_niche %in% c("Scavenger", "Vertivore"), "carnivore",
        trophic_niche %in% c("Aquatic predator",  "Herbivore aquatic"), "aquaticFeeder",
        trophic_niche %in% c("Frugivore", "Granivore", "Nectarivore"), "cropFeeder",
        default = "terrestrialFeeder"),
      altitude_cat = fcase(
        minimum_altitude < 500, "lowland",
        minimum_altitude >= 500 & minimum_altitude < 1500, "midland",
        minimum_altitude >= 1500, "highland"
      )
    )][, ..traits_select]
  
  wnv_dt <- merge(wnv_dt, traits_dt, by = "birdlife_name")
  
  
  # keep only the data that is complete
  wnv_dt[
    , data_complete := apply(.SD, 1, function(row) all(!is.na(row) & row != ""))]
  wnv_dt <- wnv_dt[data_complete == T][, data_complete := NULL]
  
  # and minimum 20 tested individuals per species
  wnv_dt[, ind_per_sp := sum(total_tested), by = avilist_name]
  wnv_dt <- wnv_dt[ind_per_sp >= 20]
  
  big_tests <- c("serological_ELISA", "molecular_detection", "serological_VNT")
  wnv_dt[, method_cat := fifelse(!method %in% big_tests, "other", method)]
  
  wnv_dt[, ':=' (
    freshwater = relevel(factor(freshwater), ref = "other"),
    migration = relevel(factor(migration), ref = "sedentary"),
    primary_lifestyle = relevel(factor(primary_lifestyle), ref = "other"),
    trophic_niche = relevel(factor(trophic_niche), ref = "terrestrialFeeder"),
    altitude_cat = relevel(factor(altitude_cat), ref = "lowland"),
    habitat = relevel(factor(habitat), ref = "greenLandscape"),
    sociality = relevel(factor(sociality), ref = "nonSocial")
  )]
  
  # save data
  fwrite(wnv_dt, file.path(cluster_dir, "1_WNV_prevalence_data_model.csv"))
  
  #PHYLOGENETIC AUTOCORRELATION
  bird_tree <- readRDS(file.path(data_dir, "00_bird_tree_match_avilist.rds"))
  bird_tree <- drop.tip(
    bird_tree,
    bird_tree$tip.label[!bird_tree$tip.label %in% unique(wnv_dt$avilist_name)]
    )
  # this way we generate covariance matrix that will be inputted in the model as:
  # gr(scientific_name_bt, cov = A) - we make sure that species are correlated
  # as specified by the covariance matrix A
  A <- vcv.phylo(bird_tree)
  saveRDS(A, file.path(cluster_dir, "1_bird_tree_A.rds"))
  
}

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


# 1: check distributions -----------------------------------------------------

# checking which variables are correlated
# log_var <- grep("^log_|_log$", names(wnv_dt), value = T)
# 
# psych::pairs.panels(wnv_dt[, ..log_var],
#              method = "pearson",
#              hist.col = "#00AFBB",
#              density = TRUE)

full_f <- "positive | trials(total_tested) ~ 1 + tarsus_log + longevity_log + clutch_max_log + log_human_population_density + abundance_log + freshwater + migration + sociality + primary_lifestyle + nest_placement + trophic_niche + altitude_cat + (1 | gr(avilist_name,cov=A)) + (1 | method_cat) + (1 | country:sampling_year)"


# get the log warning files
logfile <- file.path(model_dir, "1_check_distributions_warnings.log")
logcon <- file(logfile, open = "wt")

distributions <- list(
  binom = binomial(), 
  beta_binom = beta_binomial(), 
  zero_binom = zero_inflated_binomial(), 
  zero_beta_binom = zero_inflated_beta_binomial()
)

for(dn in names(distributions)){

    # TRY
    mname <- paste0("1_", dn, "_full_model.rds")
    
    # convert to brms formula object
    f <- bf(formula(full_f), family = distributions[[dn]])

    # Fit the model only if not already saved
    if (!file.exists(mname)) {
      
      withCallingHandlers({
        m <- brm(
          data = wnv_dt,
          data2 = list(A = A),
          formula = f,
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

  
