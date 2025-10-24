library(data.table)
library(psych)
library(brms)
library(vegan)

cluster_dir <- here::here("Data_for_cluster")

# 0: Load data cluster -------------------------------------------------------

A <- readRDS(file.path(cluster_dir, "1_bird_tree_A.rds"))

wnv_dt <- fread(file.path(cluster_dir, "1_WNV_prevalence_data_model.csv"))

wnv_dt[, prop := positive/total_tested*100]

# selected traits
traits_select <- c(
  "birdlife_name", "mass_log", "tarsus_log", "hwi_log", "longevity_log",
  "clutch_max_log", "log_human_population_density", "abundance_log",
  "habitat", "habitat_TP2019", "freshwater", "migration", "sociality",
  "primary_lifestyle", "nest_placement", "trophic_niche", "altitude_cat"
)


#first, lets explore the database a bit
numeric_fixed <- grep(".*_log|^log_", traits_select, value = T)

pairs.panels(wnv_dt[, ..numeric_fixed],
             method = "pearson",
             hist.col = "#00AFBB",
             density = TRUE)


#let's do a PCA to see how these variables cluster
wnv_dt[, coord_dec := paste(long, lat, sep = ", ")]
loc_group <- wnv_dt[, c("coord_dec", numeric_fixed), with = FALSE]
loc_group_mean <- loc_group[
  , lapply(.SD, mean, na.rm = TRUE), by = coord_dec, .SDcols = numeric_fixed]

pca_ord <- rda(loc_group_mean[, ..numeric_fixed], scale=TRUE)
summary(pca_ord)
b=biplot(pca_ord, xlab="PC1  41%", ylab="PC2 20%")
text(b$species[,1],b$species[,2],rownames(b$species),col="red")


