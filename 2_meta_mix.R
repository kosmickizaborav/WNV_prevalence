library(tidyverse)
library(data.table)
library(metafor)
library(ape)
library(ggtree)
library(aplot)
library(ggplot2)
source("0_helper_functions.R")
source("color_palette.R")

data_dir <- here::here("Data")
graphs_dir <- file.path(data_dir, "Graphs")

dt <- readxl::read_xlsx(
  file.path(data_dir, "WNV_Host_Susceptability_Review.xlsx"), 
  sheet = "observations") |> 
  janitor::clean_names()


# the ones that don't have bird_life_name are all at the genus level plus
# Ectopistes_migratorius that is extinct
dt <- setDT(dt)[class == "Aves"][bird_life_name != "NA"]
dt <- dt[ , .(
  species_tree = species_name, bl_name = bird_life_name, total_tested, 
  positive_m1 = positive_individuals_m1, 
  positive_m2 = positive_individuals_m2, method_1, method_2, sampling_year, 
  country, region = province_region, ref_id = ref_id_2, long, lat)]

# IF THERE IS A MORE SPECIFIC TEST THAN ELISA TAKE THOSE RESULTS
dt[, positive_m3 := fifelse(
  method_1 == "serological.ELISA" & method_2 != "NA", 
  as.numeric(positive_m2), positive_m1)][
    , method_3 := fifelse(
      method_1 == "serological.ELISA" & method_2 != "NA", 
      method_2, method_1
    )]


# rename to new birdlife and add phylogeny
# there was 10 entries that didn't match with any birdlife_name and all "Corvus"
dt <- rename_to_birdlife(dt, species_name = "bl_name")[!is.na(birdlife_name)][
  , name_type := NULL]
dt <- add_birdlife_phylogeny(dt, species_name = "birdlife_name")



dt[, country_positive := sum(positive_m3) > 0, by = country]
dt[, ref_n := uniqueN(ref_id), by = birdlife_name]

bird_db <- dt[total_tested >= 5 & country_positive]


wnv_summary <- bird_db |> 
  mutate(
    m1 = positive_m1/total_tested, 
    weights = 1 + log10(total_tested)
  ) |> 
  summarize(
    n_studies = length(unique(ref_id)), 
    positive_ind_m1 = sum(positive_m3),
    total_ind = sum(total_tested), 
    mean_m1 = mean(m1),
    w_mean_m1 = weighted.mean(m1, w = weights),
    pooled_m1 = positive_ind_m1/total_ind,
    #pooled_m1_w = pooled_m1 * (1 + log10(total_ind)),
    family = unique(family),
    order = unique(order),
    .by = birdlife_name
  ) 


# calculating effect size for each study
# the double-arcsine transformation because there is many 0 and 1 
# targ can also be set to list(ni = 1/(pes.da$se)^2) 
meta <- escalc(
  xi = positive_m3, 
  ni = total_tested, 
  measure = "PFT", 
  add = 0, 
  data = bird_db
)


# Assumption 1: Do not assume a common between-study variance component 
# (do not pool within-group estimates of between-study variance)
# we first fit a random-effects model for each subgroup and
# then we combine the estimated statistics into a data frame

# fitting model for each species
pes.da_all <- meta |> 
  group_split(birdlife_name) |> 
  map(~ {
    
    list(
      pes.da = rma(
        yi, vi, 
        data = .x, 
        method = "REML", 
        control = list(stepadj = 0.5), # added to converge
        level = 95  # level = 99 because the value of tau i i2 falls out of CI
        # reml better for the small sample size
        # "REML" = the restricted maximum-likelihood estimator
        # alternative and default "DL"
      ),
      
      birdlife_name = unique(.x$birdlife_name)
      
    )
    
  }, 
  .progress = T
  )

pes_all <- pes.da_all |> 
  map(~ {
    
    pes.da <- .x[["pes.da"]]
    
    pes <- predict(
      pes.da, 
      transf = transf.ipft.hm,
      targ = list(ni=1/(pes.da$se)^2)
      # alternative targ list(ni = dat$total) for harmonic mean
    )
    
    tibble(
      birdlife_name = .x[["birdlife_name"]],
      pred = pes$pred, 
      ci.lb = pes$ci.lb, 
      ci.ub = pes$ci.ub, 
      pi.lb = pes$pi.lb, 
      pi.ub = pes$pi.ub
    )
    
    # using 1/̄ v as the total sample size estimate
    # Instead of the harmonic mean, Barendregt et al. (2013) and 
    # Xu et al. (2021) recommend using 1/v (variance) as the 
    # estimate for the total sample size, because harmonic mean
    # becomes numerically unstable when sin t is close to 0 or 1,
    # leading to potentially misleading results. (p. 14)
    
  }
  ) |> 
  bind_rows()

# dat.diffvar <- pes.da_all |> 
#   map(~ {
#     
#     pes.da <- .x[["pes.da"]]
#     
#     tibble(
#       estimate = pes.da$b,
#       stderror = pes.da$se, 
#       moderator = unique(.x[["scientific_name"]]), 
#       tau2 = round(pes.da$tau2, 4)
#     )
#     
#   }) |> 
#   bind_rows()
# 
# pes.da_no_subgroups <- rma(
#   yi, 
#   vi, 
#   data = meta, 
#   method = "REML", 
#   mods = ~ factor(scientific_name) - 1
#   )
# 
# 
# pes_no_subgroups <- predict(
#   pes.da_no_subgroups, 
#   transf = transf.ipft.hm,
#   targ = list(ni=1/(pes.da_no_subgroups$se)^2)
#   )
# 
# 
# pes.moderator_same_var <- rma(
#   estimate, 
#   sei = stderror, 
#   method = "REML", 
#   data = dat.diffvar
#   )
# 
# 
# effect_glmm <- metaprop(
#   data = bird_db, 
#   event = positive_study_species, 
#   n = total_study_species, 
#   method = "GLMM", # alternative method = "Inverse
#   #incr = 0, #important for inverse 
#   studylab = bird_db$author_year, 
#   random = TRUE,
#   sm = "PLOGIT", # double arcsine sm = "PFT", 
#   # backtransf = TRUE, 
#   subgroup = factor(scientific_name), 
#   control = list(stepadj = 0.5)
# )


# # second step of subgroup meta-analysis fixed effects ---------------------
# 
# # second step
# # we fit a fixed-effect model to compare the two estimated logit transformed 
# # proportions and recalculate the summary proportion.
# subganal.moderator <- rma(
#   estimate,
#   sei = stderror,
#   mods = ~ moderator,
#   method = "FE",
#   data = dat.diffvar
#   )
# # 
# 
# pes.da.moderator <- rma(
#   estimate,
#   sei = stderror,
#   method = "FE",
#   data = dat.diffvar
#   )
# 
# pes.moderator <- predict(
#   pes.da.moderator, 
#   transf = transf.ipft.hm, 
#   targ = list(ni=1/(pes.da.moderator$se)^2) 
#   )
# 
# # display summary effect size for all subgroups (199 species)
# pes_all |> 
#   map(~{ print(.x) })
# # display subgroup analysis results print(pes.moderator)
# print(subganal.moderator) 
# #display recomputed summary effect size
# print(pes.moderator) 



# plotting graph for subgroups --------------------------------------------

all_measures <- pes_all |> 
  left_join(wnv_summary) 

fwrite(setDT(all_measures), file.path(data_dir, "2_wnv_all_measures.csv"))

all_measures <- pes_all |> 
  left_join(wnv_summary) |> 
  pivot_longer(
    cols = all_of(c("pred", "mean_m1", "w_mean_m1",  "pooled_m1")), #"ci.lb", "ci.ub",
    names_to = "variable",
    values_to = "value"
  ) |>
  mutate(
    variable = case_when(
      variable == "pred" ~ "species meta-analysis",
      variable == "mean_m1" ~ "species mean",
      variable == "w_mean_m1" ~ "weighted species mean",
      variable == "pooled_m1" ~ "pooled proportion",
      variable == "ci.lb" ~ "95% CI - lower limit",
      variable == "ci.ub" ~ "95% CI - upper limit",
    )
  )


# species level phylogeny according to BirdTree taxonomy
phy <- read.nexus(file.path(
  data_dir, "Traits_data", "HackettStage1_0001_1000_MCCTreeTargetHeights.nex"))
phy$tip.label <- gsub("_", " ", phy$tip.label)


dropTips <- setdiff(phy$tip.label, all_measures$birdlife_name)
phy_selected <- drop.tip(phy,dropTips)


p <- ggtree(phy_selected) 
p1 <- p %<+% all_measures[, c("birdlife_name", "order")] + 
  # geom_tiplab(aes(label = paste0("italic('", label, "')")), 
  #             align = TRUE, size = 5) +
  geom_tippoint(aes(colour = order)) +
  scale_color_manual(values = ord_col) +
  theme(legend.position = "none")


p2 <- all_measures |> 
  filter(birdlife_name %in% phy$tip.label) |> 
  ggplot() +
  geom_errorbar(
    aes(y = birdlife_name, 
        xmin = ci.lb, xmax = ci.ub), linewidth = 0.1
  ) +
  geom_point(
    aes(y = birdlife_name, x = value, color = variable, shape = variable)
  ) +
  geom_text(
    aes(y = birdlife_name, x = 1, label = str_c(total_ind, "|", n_studies)), size = 1.5
  ) +
  scale_shape_manual(values = c(15, 1, 3, 16)) + 
  scale_color_manual(values = c("royalblue", "orange3", "black", "red2")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "prevalence", y = "species") +
  guides(color = guide_legend(nrow = 4), shape = guide_legend(nrow = 4))


options("aplot_guides" = "keep")


p2  %>% insert_left(p1) 


ggsave(file.path(graphs_dir, "1_prev_methods.png"), height = 20, width = 15)




