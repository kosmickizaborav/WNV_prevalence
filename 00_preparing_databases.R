library(data.table)
library(ggplot2)
library(clootl)
source("0_helper_functions.R")

data_dir <- here::here("Data")
tdir <- file.path(data_dir, "Traits_data")


# WNV database ------------------------------------------------------------

dt <- readxl::read_xlsx(
  file.path(data_dir, "WNV_Host_Susceptability_Review.xlsx"), 
  sheet = "observations") |> 
  janitor::clean_names()

# the ones that don't have bird_life_name are all at the genus level plus
# Ectopistes_migratorius that is extinct
dt <- setDT(dt)[class == "Aves"][bird_life_name != "NA"]
dt <- dt[ , .(
  species_wnv = species_name, bird_life_name, e_bird_name, avibase_name, 
  total_tested, 
  positive_m1 = positive_individuals_m1, 
  positive_m2 = positive_individuals_m2, method_1, method_2, sampling_year, 
  country, region = province_region, ref_id = ref_id_2, long, lat)]

dt <- rename_to_birdlife(dt, species_name = "bird_life_name")

dt <- dt[, name_type := NULL]

# check if birdlife species is present in avilist
dt_birdlife <- match_to_avilist(
  dt[, .(bird_life_name)], species_name = "bird_life_name", add_phylo = T)

# if no birdlife match was found check ebird
dt_ebird <- dt[
  bird_life_name %in% dt_birdlife[is.na(avilist_name), bird_life_name], 
  .(e_bird_name, bird_life_name)][!is.na(e_bird_name)]
# some ebird names had specified subspecies
dt_ebird[
  , e_bird_name := sub("^([A-Za-z]+ [a-z]+).*", "\\1", e_bird_name)]
dt_ebird <- match_to_avilist(
  dt_ebird, species_name = "e_bird_name", add_phylo = T)

dt_birdlife <- unique(dt_birdlife[!is.na(avilist_name)])
dt_birdlife[, name_type := paste0("birdlife_to_", name_type)]

dt_ebird <- unique(dt_ebird[!is.na(avilist_name)])
dt_ebird[, name_type := paste0("ebird_to_", name_type)]

# last option, check avibase name
dt_avibase <- dt[
  !bird_life_name %in% c(dt_ebird$bird_life_name, dt_birdlife$bird_life_name), 
  .(bird_life_name, avibase_name)][!is.na(avibase_name)]

dt_avibase <- match_to_avilist(
  dt_avibase, species_name = "avibase_name", add_phylo = T)
dt_avibase <- unique(dt_avibase[!is.na(avilist_name)])
dt_avibase[, name_type := paste0("avibase_to_", name_type)]

avilist_dt <- rbindlist(list(dt_birdlife, dt_ebird, dt_avibase), fill = T)
avilist_dt <- avilist_dt[
  , .(bird_life_name, avilist_name, avilist_family = family, 
      avilist_order = order, name_type)]

dt <- merge(
  dt, avilist_dt, 
  by = "bird_life_name", 
  all.x = TRUE
  )

# matches to send to Alex
alex_match <- unique(dt[
  , .(avilist_name, name_type, bird_life_name, e_bird_name, avibase_name)])
fwrite(alex_match, file.path(data_dir, "rename_to_avilist.csv"))

# the only entries that don't have avilist match were "Corcus" without species
dt <- dt[!is.na(avilist_name)]

# IF THERE IS A MORE SPECIFIC TEST THAN ELISA TAKE THOSE RESULTS
dt[, positive_m3 := fifelse(
  method_1 == "serological.ELISA" & method_2 != "NA", 
  as.numeric(positive_m2), positive_m1)][
    , method_3 := fifelse(
      method_1 == "serological.ELISA" & method_2 != "NA", 
      method_2, method_1
  )]

dt[, method_3 := gsub("_|\\.", " ", method_3)]

dt <- dt[, .(
  species_wnv, avilist_name, avilist_family, avilist_order, 
  birdlife_name, name_type, 
  positive = positive_m3, 
  total_tested, 
  method = method_3, 
  sampling_year, country, region, ref_id, long, lat
)]



fwrite(dt, file.path(data_dir, "00_WNV_prevalence_data.csv"))

rm(alex_match, dt_avibase, dt_birdlife, dt_ebird)

target_sp <- dt[, unique(birdlife_name)]

# species tracking
dt_species <- unique(dt[, .(species_wnv, avilist_name, birdlife_name)])


# BIRDBASE ----------------------------------------------------------------

# rr = restricted range: Restricted range (global range size<50,000 km2) >>> 1
# lat - latitudinal range: 1 - tropical, 2 - tropical_temperate, 3 - temperate, 
#    4 - temperate-polar, 5 tropical-polar
# average mass = the average value across males, females, 
#    and unsexed individuals in the dataset
# NormMin = typical elevation limit of each species' center of abundance, 
#    normal lowest elveation of the species (L	Lowland. 0-500 m a.s.l, 
#    F	Foothill. 501-1000 m a.s.l. M	Montane. >>1000 m)
# main habitat types
# primary habitat
# primary diet 
# SOCIAL BEHAVIOUR: 
# social 1: colonial
# social 2: social (large number of birds, mixed species flocks, seasonal flocks)
# social 3: pairs (breeding pairs or small family groups)
# social 4: singly and pairs
# social 5: solitary
# social 6: lekking
# Coop: COOPERATIVE BREEDING: 0 - solitary breeder, 2 - family member helper, 
#   3 occasional help
# Brd1, Brd2 - The minimum and maximum (respectively) number of broods 
#   for the species.
# Clutch_Min, Max - Lower and upper number of eggs laid (excluding rare 
#   clutch sizes)
# Prod 1, Prod2 - Yearly productivity per pair (# young raised to fledging)
# Flightlessness - yes, no, partial
# Mig - Migration. 1 - Yearly, regular, seasonal (larger scale), and/or 
#   long-distance movement. 2 indicates partial migrant. 
# Alt - Altitudinal/elevational migration. Regular movements from high to low 
#   altitudes or vice versa based on seasons or post-breeding dispersal. 
#   Includes notable altitudinal movements on a daily basis for foraging 
#   (i.e. Barred Parakeet, Tepui Parrotlet)
# Disp - Long-distance dispersal. Usually once in a lifetime, 
#   after leaving nest (i.e. post-breeding or post-fledging dispersal). 
# Sed - sedentary


cols_of_interest <- c(
  birdlife_name = "hbw_bird_life_international_v9_1", 
  range_restriction = "rr", latitudinal_range = "lat",
  "average_mass", elev_norm_min = "norm_min", 
  "primary_habitat", "primary_diet", 
  colonial = "social_1", social = "social_2", pairs = "social_3", 
  singly = "social_4", solitary = "social_5", lekking = "social_6", 
  cooperative_breeding = "coop", max_n_broods = "brd2", "clutch_max",
  yearly_productivity_max = "prod2", "flightlessness", migration_birdbase = "mig",
  altitudinal_migration = "alt", long_distance_dispersal = "disp",
  sedentary = "sed"
)

birdbase <- readxl::read_xlsx(
  file.path(tdir, "BIRDBASE v2025.1 Sekercioglu et al. Final.xlsx"), 
  skip = 1,
  sheet = "Data") |> 
  janitor::clean_names() |> 
  dplyr::select(all_of(cols_of_interest))

setDT(birdbase)

birdbase <- birdbase[birdlife_name %in% target_sp]

# AVONET ----------------------------------------------------------------

# all 1426 species matched - all species

# all variables 
# beak: length_culmen, length_nares", width, depth, 
# tarsus: length,
# wing: length, kipps distance, secondary, hwi,
# tail: lengh
# mass:
# habitat: density
# migration
# trophic: level, niche, primary lifestyle
# range: max lat/lon, size, centroid

trt_cols <- c(
  "species1", "hand_wing_index", "mass", "tarsus_length", 
  "habitat", "habitat_density", "migration", "trophic_niche", 
  "primary_lifestyle", "range_size")

avonet <- fread(file.path(tdir, "AVONET", "AVONET1_BirdLife.csv")) |> 
  janitor::clean_names() 

# take only traits of interest
avonet <- avonet[, ..trt_cols]

avonet <- rename_to_birdlife(avonet, "species1")

avonet <- avonet[birdlife_name %in% target_sp][
  , name_type := NULL]

dt_species <- merge(
  dt_species, 
  avonet[, .(birdlife_name, avonet_sp = species1)], 
  by = "birdlife_name",
  all.x = TRUE)

avonet[, species1 := NULL]


# Tobias & Pigot 2019 -----------------------------------------------------

# 1382/1426 species matched
# missing 32
# all variables: "redlist_category", "threat", "log_range_size", 
# "log_body_mass", "diet", "foraging", "migration", "mating_system", 
# "nest_placement", "territoriality", "habitat", "island_dwelling", 
# "log_clutch_size", "log_night_lights", "log_human_population_density"

# columns of interest
trt_cols <- c(
  "species", "nest_placement", "territoriality", "mating_system", 
  "log_clutch_size", "foraging", "habitat", "log_human_population_density"
  )

tobias <- fread(file.path(tdir, "Tobias&Pigot2019.csv")) |> 
  janitor::clean_names() 

tobias <- tobias[, ..trt_cols]

setnames(tobias, old = "habitat", new = "habitat_TP2019")

tobias[, species := gsub("_", " ", species)]

tobias <- rename_to_birdlife(tobias, species_name = "species")[
  , name_type := NULL]

tobias <- tobias[birdlife_name %in% target_sp]

# removing the duplicates that were caused by synonyms in birdlife list
dupl_sps <- tobias[duplicated(birdlife_name), unique(birdlife_name)]
tobias[, duplicates := birdlife_name %in% dupl_sps]
  
tobias <- tobias[
  duplicates == F | 
    duplicates == T & birdlife_name == species | species == "Sylvia curruca"]

dt_species <- merge(
  dt_species, 
  tobias[, .(birdlife_name, tobias2019_sp = species)], 
  by = "birdlife_name",
  all.x = TRUE)

tobias[, species := NULL][, duplicates := NULL]


# Bird 2020 ---------------------------------------------------------------

# 1426/1426 found - all

col_bird4 <- c("scientific_name", "adult_survival", "age_at_first_breeding", 
               "maximum_longevity", "gen_length")

col_bird3 <- c("scientific_name", "mean_clutch_size", "minimum_altitude", 
               "maximum_altitude", "nocturnal", "freshwater")

bird3 <- fread(file.path(tdir, "Bird2020_table3.csv")) |> 
  janitor::clean_names()

bird3 <- bird3[, ..col_bird3]

bird3 <- rename_to_birdlife(bird3, species_name = "scientific_name")[
  , name_type := NULL]

# there were two entries for Ardea pacifica in bird3
# one with typo (Arden pacifica), 
# but different values for the clutch size and altitude
# i kept the one that was Ardea pacifica


bird3 <- bird3[!is.na(birdlife_name)][scientific_name != "Arden pacifica"]

dt_species <- merge(
  dt_species, 
  bird3[, .(birdlife_name, bird2020_table3_sp = scientific_name)], 
  by = "birdlife_name",
  all.x = TRUE)

bird3[, scientific_name := NULL]

bird <- fread(file.path(tdir, "Bird2020_table4.csv")) |> 
  janitor::clean_names()

bird <- bird[, ..col_bird4]
bird <- rename_to_birdlife(bird, species_name = "scientific_name")[
  , name_type := NULL]

dt_species <- merge(
  dt_species, 
  bird[, .(birdlife_name, bird2020_table4_sp = scientific_name)], 
  by = "birdlife_name",
  all.x = TRUE)

bird[, scientific_name := NULL]
bird <- bird[!is.na(birdlife_name)]

bird <- merge(bird3, bird, by = "birdlife_name", all = TRUE)
rm(bird3)

bird <- bird[birdlife_name %in% target_sp]


# Callaghan et al. 2021 ---------------------------------------------------

col_trt <- c("scientific_name", "abundance_estimate")

call2021 <- fread(file.path(tdir, "Callaghan2021.csv")) |> 
  janitor::clean_names()

call2021 <- call2021[, ..col_trt]
call2021[, scientific_name := gsub("_", " ", scientific_name)]

call2021 <- rename_to_birdlife(call2021, species_name = "scientific_name")[
  , name_type := NULL]
 
call2021 <- call2021[birdlife_name %in% target_sp]

# removing the duplicates that were caused by synonyms in birdlife list
dupl_sps <- call2021[duplicated(birdlife_name), unique(birdlife_name)]
call2021[, duplicates := birdlife_name %in% dupl_sps]

call2021 <- call2021[
  duplicates == F | 
    duplicates == T & birdlife_name == scientific_name | 
    scientific_name == "Sylvia curruca"]

dt_species <- merge(
  dt_species, 
  call2021[, .(birdlife_name, call2021_sp = scientific_name)], 
  by = "birdlife_name",
  all.x = TRUE)

call2021[, scientific_name := NULL][, duplicates := NULL]


# Sheard 2023 nest structure ----------------------------------------------


trt_cols <- c("species_scientific_name", "nest_height_min_m",
              "nest_height_max_m", "source")

# 1180/1426 species 

nest_height <- fread(
  file.path(tdir, "Sheard2023_nest_structure_Dataset-S1.csv")) |> 
  janitor::clean_names() 

nest_height <- nest_height[, ..trt_cols]

# removing no_info fields and converting them to NA
nest_height <- melt(
  nest_height,
  id.vars = c("species_scientific_name", "source"),
  measure.vars = c("nest_height_min_m", "nest_height_max_m"),
  variable.name = "nest_height"
)

# remove NA and no info, after NA will reappear if there is some data for the 
# species ay least
nest_height <- nest_height[!(is.na(value) | value == "no info")]
# back to wide format
nest_height <- dcast(
  nest_height, 
  species_scientific_name + source ~ nest_height, 
  value.var = "value"
)

# check if the value was labelled as uncertain
nest_height[, ":=" (
  uncertain_min = grepl("Uncertain", nest_height_min_m, ignore.case = T), 
  uncertain_max = grepl("Uncertain", nest_height_max_m, ignore.case = T)
)]
# if there were uncertainty in any of the measurements label the species
nest_height[, uncertainty_nest_height := fcase(
    sum(uncertain_min) >= 1 & sum(uncertain_max) >= 1, "both", 
    sum(uncertain_min) >= 1, "min",
    sum(uncertain_max) >= 1,  "max", 
    default = "none"
  ),
  by = species_scientific_name]

# for both min and max height, 
# remove uncertainty and just take the value that is given
nest_height[, (trt_cols[2:3]) := lapply(
  .SD, function(x) gsub(".*\\((\\d+)\\).*|^(\\d+)$", "\\1\\2", x)), 
  .SDcols = trt_cols[2:3]]

nest_height <- nest_height[, .(
  nest_height_min_m = min(nest_height_min_m, na.rm = T), 
  nest_height_max_m = max(nest_height_max_m, na.rm = T), 
  uncertainty_nest_height = unique(uncertainty_nest_height)), 
  by = species_scientific_name]

nest_height <- rename_to_birdlife(
  nest_height, species_name = "species_scientific_name")[, name_type := NULL]

dt_species <- merge(
  dt_species, 
  nest_height[, .(birdlife_name, sheard2023_sp = species_scientific_name)], 
  by = "birdlife_name",
  all.x = TRUE)

nest_height[, species_scientific_name := NULL]

nest_height <- nest_height[birdlife_name %in% target_sp]

# 1361/1426 found
# 
# nest_location <- fread(
#   file.path(tdir, "Sheard2023_nest_structure_Dataset-S2.csv")) |> 
#   janitor::clean_names() 
# 
# trt_cols <- c(
#   "species_scientific_name", grep("^loc", names(nest_location), value = T))
# 
# nest_location <- nest_location[, ..trt_cols]
# 
# nest_location <- rename_to_birdlife(
#   nest_location, species_name = "species_scientific_name")[
#     , species_scientific_name := NULL][, name_type := NULL]
# 
# nest_location <- nest_location[birdlife_name %in% target_sp]
# 
# # remove all rows that have NA value across all columns
# nest_location <- nest_location[
#   !apply(is.na(nest_location[, -1, with = FALSE]), 1, all)]
  

# full_dataset_traits  ----------------------------------------------------

# Create a list of your data.tables
dt_list <- list(avonet, tobias, bird, call2021, nest_height, birdbase)

# Use Reduce to iteratively left join by "birdlife_name"
final <- Reduce(
  function(x, y) merge(x, y, by = "birdlife_name", all.x = TRUE), dt_list)


rm(birdbase, dt_list, avonet, tobias, bird, call2021, nest_height, nest_location)

final[, habitat_density := fcase(
  habitat_density == 1, "dense habitat", 
  habitat_density == 2, "semi-open habitat",
  habitat_density == 3, "open habitat"
)] 

final[, migration := fcase(
  migration == 1, "sedentary", 
  migration == 2, "partially migratory", 
  migration == 3, "migratory"
)]


fwrite(final, file.path(data_dir, "00_traits_filtered.csv"))


# Selected traits --------------------------------------------------------------

birdbase <- c(
  "range_restriction", "latitudinal_range", "average_mass", "norm_min", 
  "primary_habitat", "primary_diet", "colonial", "social", "pairs", "singly", 
  "solitary", "lekking", "cooperative_breeding", "max_n_broods", "clutch_max",
  "yearly_productivity_max", "flightlessness", "migration_birdbase", 
  "altitudinal_migration", "long_distance_dispersal", "sedentary"
)

# AVONET
avonet <- c(
  "hand_wing_index", "mass", "tarsus_length", 
  "habitat", "migration", "trophic_niche", "primary_lifestyle", "range_size")

# Tobias & Pigot 2019
tobias_pigot_2019 <- c(
  "species", "nest_placement", "territoriality", "mating_system", 
  "log_clutch_size", "foraging", "habitat_tobias2019", 
  "log_human_population_density"
)

# bird2020
bird2020 <- c("adult_survival", "age_at_first_breeding", 
              "maximum_longevity", "gen_length", 
              "mean_clutch_size", "minimum_altitude", 
              "maximum_altitude", "nocturnal")

callaghan2021 <- "abundance_estimate"

sheard2023 <- c(
  "nest_height_min_m", "nest_height_max_m", "uncertainty_nest_height"
)


traits_table <- data.table(
  trait = c(
    avonet, 
    tobias_pigot_2019, 
    bird2020, 
    callaghan2021, 
    sheard2023
  ),
  source = c(
        rep("BIRDBASE", length(names(birdbase)),
        rep("AVONET", length(avonet)),
        rep("Tobias & Pigot, 2019", length(tobias_pigot_2019)),
        rep("Bird, 2020", length(bird2020)),
        rep("Callaghan, 2021", length(callaghan2021)),
        rep("Sheard, 2023", length(sheard2023))
    )
  )
)
  
  fwrite(traits_table, file.path(data_dir, "00_traits_source-variable.csv"))
  
  

# PLOT: data availability ------------------------------------------------------

summary <- melt(
  rbind(
    final[, lapply(.SD, function(x) sum(!is.na(x)))][, type := "present"],
    final[, lapply(.SD, function(x) sum(is.na(x)))][, type := "missing"]
  ),
  id.vars = "type",
  variable.name = "name",
  value.name = "value"
)

summary[, name := gsub("_", " ", name)]

brks <- seq(1, nrow(final), 100)

summary |> 
  ggplot() +
  geom_col(aes(x = name, y = value, fill = type)) +
  geom_hline(yintercept = brks, color = "gray33", linetype = "dashed") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), breaks = brks) +
  labs(
    title = "Availability of trait data", 
    y = "Number of species", 
    x = 'Trait', 
    fill = ""
  ) +
  theme(legend.position = "bottom") + 
  coord_flip()

ggsave(here::here("Data", "Graphs", "00_trait_data_availability.png"),
  width = 8, height = 6, units = "in"
)



# Bird tree data -----------------------------------------------------------

full_tree <- extractTree()

sp_tree <- data.table(
  species_tree = gsub("_", " ", full_tree$tip.label))

sp_tree <- match_to_avilist(sp_tree, species_name = "species_tree")

sp_tree <- sp_tree[avilist_name %in% avilist_dt$avilist_name]

pruned_tree <- extractTree(species = unique(sp_tree$species_tree))

dt_species <- merge(
  dt_species, 
  sp_tree[, .(avilist_name, species_tree)], 
  by = "avilist_name",
  all.x = TRUE)

# rename tree labels to AVILIST names
sp_tree <- sp_tree[
  , names(sp_tree) := lapply(.SD, function(x) gsub(" ", "_", x))]

pruned_tree$tip.label <- sp_tree[
  match(pruned_tree$tip.label, species_tree), avilist_name]

saveRDS(pruned_tree, file.path(data_dir, "00_bird_tree_match_avilist.rds"))


# crosswalk scientific names all databases --------------------------------

fwrite(
  dt_species, file.path(data_dir, "00_species_crosswalk_all_databases.csv"))
