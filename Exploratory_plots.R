library(data.table)
library(ggplot2)
library(colorspace)
library(MetBrewer)
library(patchwork)
library(ggtree)
library(ggtreeExtra)
library(ape)
library(sf)
source("0_helper_functions.R")
source("color_palette.R")

data_dir <- here::here("Data")
graphs_dir <- file.path(data_dir, "Graphs")
dir.create(graphs_dir, showWarnings = F)

traits <- fread(file.path(data_dir, "00_traits_filtered.csv")) |> 
  janitor::clean_names()

dt <- readxl::read_xlsx(
  file.path(data_dir, "WNV_Host_Susceptability_Review.xlsx"), 
  sheet = "observations") |> 
  janitor::clean_names()

sps <- unique(dt$species_name)

# the ones that don't have bird_life_name are all at the genus level plus
# Ectopistes_migratorius that is extinct
dt <- setDT(dt)[class == "Aves"][bird_life_name != "NA"]
dt <- dt[ , .(
  species_tree = species_name, bl_name = bird_life_name, e_bird_name, total_tested, 
  positive_m1 = positive_individuals_m1, 
  positive_m2 = positive_individuals_m2, method_1, method_2, sampling_year, 
  country, region = province_region, ref_id = ref_id_2, long, lat)]

dt <- match_to_avilist(dt, species_name = "bird_life")

check <- dt[is.na(avilist_name)][, c("avilist_name", "name_type") := NULL]

check <- match_to_avilist(check, species_name = "e_bird_name")

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


# P 1: species vs. order ----------------------------------------------------------

n_sp <- uniqueN(dt$birdlife_name)
n_fam <- uniqueN(dt$family)
n_ord <- uniqueN(dt$order)

sp_per_order <- dt[, .(n_species = uniqueN(birdlife_name)), by = .(order)]

# Summarize: number of species per family within each order
sp_per_order |> 
  ggplot(
    aes(x = reorder(order, -n_species), y = n_species, fill = order)
  ) +
  geom_bar(stat = "identity", color = "gray33", alpha = 0.4, linewidth = 0.2) +
  geom_text(aes(label = n_species), vjust = -0.5) +
  labs(
    x = "Order",
    y = "Number of species",
    title = sprintf(
      "%s bird species from %s families and %s orders", n_sp, n_fam, n_ord)
  ) +
  scale_fill_manual(values = ord_col) +
  scale_y_continuous(
    expand = c(0, 0), limits = c(0, 1.1*max(sp_per_order$n_species))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold"), 
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) 

ggsave(
  filename = file.path(graphs_dir, "1_species_per_order.png"),
  width = 10, height = 6, dpi = 300
)

rm(sp_per_order)


# P 2: species per family ------------------------------------------------------

# Count unique species per family within each order
sp_per_family <- dt[
  , .(n_species = uniqueN(birdlife_name)), by = .(order, family)][
  , ':=' (
    total_species = sum(n_species),
    n_family = uniqueN(family)), 
  by = order][total_species >= 10][
  , order_family := paste(order, family, sep = "_")][
  , order_lab := paste0(order, " [N = ", n_family, "]")]


bp <- lapply(split(sp_per_family, by = "order"), function(fm_dt){
  
  setorder(fm_dt, n_species)
  fm_dt[, fm_id := .I]
  
  brks_y <- unique(fm_dt$n_species)
  lbs_y <- brks_y
  
  if(nrow(fm_dt) == 1){
    brks_y <- c(0, brks_y)
    lbs_y <- c("", lbs_y)
  } 
  if(nrow(fm_dt) > 10){
    sq <- 1:length(brks_y)
    lbs_y <- ifelse(sq %% 3 == 0, brks_y[sq], "")
  }
  
  brks_x <- fm_dt$fm_id
  lbs_x <- fm_dt$family
  
  if(nrow(fm_dt) > 5){
    stp <- ceiling(nrow(fm_dt) / 5)
    lbs_x <- ifelse(brks_x %% stp == 0, lbs_x[brks_x], "")
  }
  
  clr <- ord_col[unique(fm_dt$order)]
  
  fm_dt |> 
    ggplot() + 
    geom_bar(
      aes(x = fm_id, y = n_species, fill = order),
      stat = "identity", color = "white", alpha = 0.8, linewidth = 0.2
    ) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = ord_col) +
    scale_x_continuous(
      breaks = brks_x, labels = lbs_x
    ) +
    scale_y_continuous(
      expand = c(0.01, 0.01), breaks = brks_y, labels = lbs_y
    ) +
    facet_wrap(~order_lab) +
    theme_bw() +
    theme(
      legend.position = "none", 
      strip.background = element_rect(fill = alpha(clr, 0.5), color = "gray33"),
      strip.text = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(0,0,0,0)
    ) +
    labs(x = "", y = "")
    
})

wrap_plots(bp, plot_spacing = units(0, "mm")) + 
  plot_annotation(
    title = "Number of families and species across orders",
    subtitle = "only orders with 10 or more species shown",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
      plot.subtitle = element_text(hjust = 0.5)
    )
  )


ggsave(
  filename = file.path(graphs_dir, "1_species_per_family_pie.png"),
  width = 13, height = 10, dpi = 300
)

rm(sp_per_family)



# P 3: species per method -------------------------------------------------

m_mtd_dt <- dt[, .(ind_tot = sum(total_tested)), by = .(method_3)]
m_mtd_dt[, method_1 := gsub("_|\\.", " ", method_3)]

m_mtd_dt |> 
  ggplot( aes(x = reorder(method_3, -ind_tot), y = ind_tot, fill = method_3)) + 
  geom_bar(
    stat = "identity", color = "gray33", alpha = 0.8, linewidth = 0.2
  ) +
  geom_text(aes(label = ind_tot), vjust = -0.5) +
  labs(
    x = "Method",
    y = "Number of individuals tested",
    title = "Number of individuals tested per method"
  ) +
  scale_fill_manual(values = MetBrewer::met.brewer("Hiroshige", n = nrow(m_mtd_dt))) +
  scale_y_continuous(
    expand = c(0, 0), limits = c(0, 1.1*max(m_mtd_dt$ind_tot))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold"), 
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  filename = file.path(graphs_dir, "1_tested_per_method3.png"),
  width = 8, height = 6, dpi = 300, bg = "white"
)

rm(m_mtd_dt)

# P 4: individuals per species -------------------------------------------------

ind_dt <- dt[, .(ind_tot = sum(total_tested)), by = .(birdlife_name, order)]
ind_dt[, ind_tot_log := log10(ind_tot)]
setorder(ind_dt, -ind_tot)
ind_dt[, sp_id := .I]

logs <- 0:floor(max(ind_dt$ind_tot_log))
logs_labs <- 10^logs

# Define the sample sizes of interest
sample_sizes <- c(max(ind_dt$ind_tot), 100, 50, 30, 10, 5, 4, 3, 2, 1)
# For each sample size, find the max sp_id
y_brks <- ind_dt[
  ind_tot %in% sample_sizes, .(max_sp_id = max(sp_id)), by = ind_tot]

ind_dt |> 
  ggplot(aes(x = sp_id, y = ind_tot_log, color = order)) + 
  geom_point(alpha = 0.8) +
  labs(
    x = "Ordinal number of species",
    y = "Number of individuals tested\n[logarithmic scale]",
    title = "Number of individuals tested per species"
  ) +
  geom_vline(
    xintercept = y_brks$max_sp_id, linetype = "dashed", 
    color = "black", alpha = 0.5, linewidth = 0.5
  ) +
  scale_color_manual(values = ord_col) +
  scale_y_continuous(
    breaks = c(logs, log10(c(2, 3, 5, 30, 50))),
    labels = c(logs_labs, 2, 3, 5, 30, 50),
    limits = c(0, NA)
  ) +
  scale_x_continuous(
    breaks = y_brks$max_sp_id
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
    strip.text = element_text(face = "bold"), 
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  filename = file.path(graphs_dir, "1_individuals_per_species.png"),
  width = 8, height = 5, dpi = 300, bg = "white"
)

rm(ind_dt, y_brks, sample_sizes)



# P 5: Phylogenetic tree -------------------------------------------------------

smpl_min <- c(1, 10, 30)

lapply(smpl_min, function(s_min){
  
  tree <- read.nexus(file.path(
    data_dir, "Traits_data", "HackettStage1_0001_1000_MCCTreeTargetHeights.nex")) 
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  
  sp_dt <- dt[
    , .(
      p_wnv = sum(positive_m3)/sum(total_tested), 
      total = sum(total_tested)), 
    by = .(birdlife_name, order)]
  # filter minimum sample size
  sp_dt <- sp_dt[total >= s_min]
  sp_dt[, p_wnv := ifelse(p_wnv > 0, "positive", "negative")]
  
  sp_dt <- unique(sp_dt[, label := birdlife_name])[
    , .(label, p_wnv, order, total)]
  
  sp_dt <- sp_dt[label %in% tree$tip.label]
  
  tree <- drop.tip(tree, setdiff(tree$tip.label, sp_dt$label))
  
  # 4. Basic tree plot
  tree_plot <- ggtree(tree, layout="circular", color = "gray88") %<+% sp_dt
  
  stl <- "colored by the test result"
  
  # 5. Join data and color **terminal branches**
  tree_plot +
    # geom_fruit(
    #   data = sp_dt,
    #   geom = geom_tile,
    #   mapping = aes(y = label, fill = order),
    #   width = 0.1,          # thickness of the ring
    #   offset = 0.2,         # distance from the tips
    #   axis.params = list()  # remove axis
    # ) +      # Color branches by WNV status
    geom_tippoint(aes(color = p_wnv), size = 1) +
    geom_tiplab(aes(color = p_wnv), size = 1, offset = 0.5) +  
    scale_color_manual(
      values = c("positive" = "#FF9966", "negative" = "#00416A"),
      na.value = "grey66"
    ) +
    labs(
      title = "Phylogenetic tree of bird species tested for WNV",
      subtitle = ifelse(
        s_min > 1,
        paste(stl, sprintf("(only species with >= %d individuals tested)", s_min)), 
        paste(stl, "(all tested species)")),  
      color = "test result"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
    
  pname <- sprintf("1_phylo_tree_wnv_status_min%s.png", s_min)
  
  ggsave(
    filename = file.path(graphs_dir, pname),
    width = 6, height = 8, dpi = 300, bg = "white"
  )
  
  sp_dt[, .(n = .N), by = p_wnv][, n_sample := s_min]
  
})



# P 6: Worldmap -----------------------------------------------------------

# world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
# world$country <- world$sovereignt
# world$region <- world$admin

world_regions <- rnaturalearth::ne_states(returnclass = "sf")

# removing invalid geometries - that created lines across maps
# fixing vald maps
#world_regions <- sf::st_make_valid(world_regions[, c("name", "geometry")])
world_regions <- world_regions[st_is_valid(world_regions), ]

# Add bounding box data
# world_regions$max_lat <- st_coordinates(st_centroid(world_regions))[,2]
# 
# # Remove "regions" above latitude 85 (the Arctic lines)
# world_regions <- world_regions[
#   world_regions$max_lat < 50 & world_regions$max_lat > -50, ]



dt_sf <- sf::st_as_sf(dt, coords = c("long", "lat"), crs = 4326, remove = F)

dt_sf <- st_join(
  dt_sf, world_regions, join = st_nearest_feature, left = T)

# no_data <- setDT(dt_sf)[is.na(name)][, name := NULL]
# no_data <- st_as_sf(no_data, coords = c("long", "lat"), crs = 4326, remove = F)
# no_data <- st_join(
#   no_data, world_regions, join = st, left = T)

dt_sf <- setDT(dt_sf)[, .(total_tested, positive_m3, name)][
  , .(
    region_tested = sum(total_tested), 
    region_positive = sum(positive_m3)),
  by = name]

dt_sf[, region_positive := fifelse(region_positive == 0, 0.1, region_positive)]

world_virus <- merge(world_regions, dt_sf, by = "name", all.x = TRUE)

ggplot() + 
  geom_sf(data = world_virus, aes(fill = region_positive)) +
  scale_fill_gradientn(
    colors = wesanderson::wes_palette("Zissou1", 100, type = "continuous"), 
    trans = "log10", na.value = "white", name = "N positive (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(
    fill = guide_colourbar(title.position="top", title.hjust = 0.5, angle = 90)
  )


ggsave(
  filename = file.path(graphs_dir, "1_worldmap_positive_per_region.png"),
  width = 10, height = 6, dpi = 300, bg = "white"
)


ggplot() + 
  geom_sf(data = world_virus, aes(fill = region_tested)) +
  scale_fill_viridis_c(
    trans = "log10", na.value = "white", name = "N tested (log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(
    fill = guide_colourbar(title.position="top", title.hjust = 0.5, angle = 90)
  )


ggsave(
  filename = file.path(graphs_dir, "1_worldmap_tested_per_region.png"),
  width = 10, height = 6, dpi = 300, bg = "white"
)

rm(world_regions, dt_sf, world_virus)


# P 7: trait availability-------------------------------------------------------

smpl_min <- c(1, 10, 30)

lapply(smpl_min, function(s_min){
  
  traits_dt <- fread(file.path(data_dir, "00_traits_filtered.csv"))
  
  traits_souce <- fread(file.path(data_dir, "00_traits_source-variable.csv"))
  
  sp_dt <- dt[, .(total_tested = sum(total_tested)), by = birdlife_name][
    total_tested >= s_min]
  
  traits_dt <- traits_dt[birdlife_name %in% unique(sp_dt$birdlife_name)]
  
  
  summary <- melt(
    rbind(
      traits_dt[, lapply(.SD, function(x) sum(!is.na(x) & x != ""))][
        , type := "present"],
      traits_dt[, lapply(.SD, function(x) sum(is.na(x) | x == ""))][
        , type := "missing"]
    ),
    id.vars = "type",
    variable.name = "trait",
    value.name = "count"
  )
  
  summary <- merge(summary, traits_souce, by = "trait", all.x = T)
  
  summary[, trait := gsub("_", " ", trait)]
  
  setorder(summary, source, count)
  summary[, trait_id := .GRP, by = trait]
  
  brks <- seq(1, nrow(traits_dt), 100)
  
  source_dt <- summary[, .(
    xmin = min(trait_id) - 0.5,
    xmax = max(trait_id) + 0.5,
    ymin = -max(brks)*0.3,
    ymax = 0, 
    ymed = -max(brks)*0.15),
    by = source][, xmed := mean(c(xmin, xmax)), by = source]
  
  source_dt[source != "Callaghan, 2021", source := gsub(",", " et al. \n", source)]
  

  
  limits_x <- c(min(source_dt$ymin), max(brks) +1)
  
  summary |> 
    ggplot() +
    geom_rect(
      data = source_dt,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
      fill = "white", color = "gray33"
    ) +
    geom_text(
      data = source_dt,
      aes(x = xmed, y = ymed, label = source), size = 3
    ) +
    geom_col(aes(x = trait_id, y = count, fill = type), alpha = 0.5) +
    scale_fill_manual(
      values = c("present" = "#4DB5C5", "missing" = "#B53389")
    ) +
    geom_hline(yintercept = brks, color = "gray33", linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(
      breaks = summary$trait_id, labels = summary$trait, expand = c(0,0)
    ) +
    scale_y_continuous(
      expand = c(0, 0), breaks = brks
    ) +
    labs(
      title = "Availability of trait data", 
      y = "Number of species", 
      x = 'Trait', 
      fill = ""
    ) +
    theme(legend.position = "bottom") + 
    coord_flip()
  
  pname <- sprintf("1_trait_data_availability_min_%d.png", s_min)
  
  ggsave(file.path(graphs_dir, pname), width = 8, height = 6, units = "in"
  )
})




# P 8: Categorical traits ------------------------------------------------------

traits_dt <- fread(file.path(data_dir, "00_traits_filtered.csv"))
# p_dt <- fread(file.path(data_dir, "2_wnv_all_measures.csv"))

p_dt <- dt[
  , n_ind_sp := sum(total_tested), by = birdlife_name][n_ind_sp >= 10]
p_dt[, ":=" (
  p_m3 = positive_m3/total_tested, 
  weight_m3 = 1 + log10(total_tested))]

p_dt <- p_dt[
  , .(
    w_mean_m3 = weighted.mean(p_m3, weight_m3, na.rm = TRUE), 
    total_ind = sum(total_tested), 
    n_studies = uniqueN(ref_id)
  ), by = .(birdlife_name, order)]

p_dt <- merge(p_dt, traits_dt, by = "birdlife_name", all.x = TRUE)

categorical_traits <- c(
  "habitat", "primary_lifestyle", "migration","habitat_density",
  "trophic_niche", "foraging", "nocturnal", "nest_placement", "territoriality"
  
) 

lapply(categorical_traits, function(cat_var){
  # Make a copy and set cat_variable column
  plot_dt <- copy(p_dt)[, cat_variable := factor(get(cat_var))][
    , .(p = w_mean_m3, cat_variable, birdlife_name)]
  
  counts_dt <- plot_dt[
    , .(n_sp = uniqueN(birdlife_name)), by = cat_variable]
  
  # Build plot and assign to p
  ggplot(plot_dt, aes(x = cat_variable, y = p, fill = cat_variable)) +
    geom_boxplot(alpha = 0.8) +
    geom_text(
      data = counts_dt, 
      aes(x = cat_variable, y = 1.1, label = n_sp), 
      vjust = 0.5, color = "gray33", size = 3
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    MetBrewer::scale_fill_met_d("Renoir") +
    labs(
      x = gsub("_", " ", cat_var), 
      y = "weighted mean prevalence"
    )
  
  # Save plot
  ggsave(
    filename = file.path(graphs_dir, sprintf("1_prev_%s_boxplot.png", cat_var)),
    width = 11, height = 6, dpi = 300, bg = "transparent"
  )
  
})

nest_locs <- c(
  "loc_artificial", "loc_earth_hole", "loc_ground", "loc_rocks",
  "loc_tree_hole", "loc_veg", "loc_water", "w_mean_m3", "birdlife_name")

plot_dt <- copy(p_dt)[, ..nest_locs]

plot_dt <- melt(
  plot_dt, 
  id.vars = c("w_mean_m3", "birdlife_name"),
  variable.name = "nest_location", 
  value.name = "presence"
)

plot_dt <- plot_dt[
  , nest_location := gsub("loc ", "", gsub("_", " ", nest_location))][
    presence == 1]

counts_dt <- plot_dt[
  , .(n_sp = uniqueN(birdlife_name)), by = nest_location]

# Build plot and assign to p
ggplot(plot_dt, aes(x = nest_location, y = w_mean_m3, fill = nest_location)) +
  geom_boxplot(alpha = 0.8) +
  geom_text(
    data = counts_dt, 
    aes(x = nest_location, y = 1.1, label = n_sp), 
    vjust = 0.5, color = "gray33", size = 3
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  MetBrewer::scale_fill_met_d("Renoir") +
  labs(
    x = "nest location", 
    y = "weighted mean prevalence"
  )

# Save plot
ggsave(
  filename = file.path(graphs_dir, "1_prev_nest_location_boxplot.png"),
  width = 11, height = 6, dpi = 300, bg = "transparent"
)



# Numerical traits --------------------------------------------------------


numerical_traits <- c(
  "hand_wing_index", "mass", "tarsus_length", "log_clutch_size", 
  "adult_survival", "age_at_first_breeding", "maximum_longevity", 
  "gen_length", "mean_clutch_size", "minimum_altitude", "maximum_altitude",
  "abundance_estimate", "nest_height_min_m", "nest_height_max_m")

lapply(numerical_traits, function(num_var){
  print(num_var)
  
  # Make a copy and set cat_variable column
  plot_dt <- copy(p_dt)[, num_variable := get(num_var)][
    , .(p = w_mean_m3, num_variable, birdlife_name, order)]
  
  # Build plot and assign to p
  ggplot(plot_dt, aes(x = num_variable, y = p, color = order)) +
    geom_point(alpha = 0.8) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = ord_col) +
    scale_x_log10() +
    labs(
      x = gsub("_", " ", paste(num_var, "[log scale]")), 
      y = "weighted mean prevalence"
    )
  
  # Save plot
  ggsave(
    filename = file.path(graphs_dir, sprintf("1_prev_%s_boxplot.png", num_var)),
    width = 11, height = 6, dpi = 300, bg = "transparent"
  )
  
})






# numbers -----------------------------------------------------------------

prev <- dt[, .(
  positive = sum(positive_m1, na.rm = TRUE), 
  total = sum(total_tested, na.rm = TRUE), 
  prev = sum(positive_m1, na.rm = TRUE) / sum(total_tested, na.rm = TRUE)), 
  by = .(bl_name, family)][total >= 10]
cols <- cols_dt[, color]
names(cols) <- cols_dt$family

copy(dt)[, positive_m1 := ifelse(positive_m1 == 0, 0.1, positive_m1)] 

copy(dt)[, positive_m1 := ifelse(positive_m3 == 0, 0.1, positive_m3)] |> 
  ggplot() +
  geom_point(
    aes(x = total_tested, y = positive_m1, color = as.factor(order)), 
    alpha = 0.5
  ) +
  theme_bw() +
  theme(legend.position = "none") + 
  scale_color_manual(values = ord_col) + 
  scale_y_log10() +
  scale_x_log10() +
  labs(
    x = "tested", 
    y = "positive",
    title = "Positive vs. Tested") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

ggsave(
  filename = file.path(graphs_dir, "1_all_positive_tested_log10.png"),
  width = 6, height = 6, dpi = 300, bg = "transparent"
)

copy(dt)[, .(
  positive = sum(positive_m1), 
  total = sum(total_tested)), 
  by = .(birdlife_name, family)] |> 
  ggplot() +
  geom_point(
    aes(x = total, y = positive, color = as.factor(family)), alpha = 0.5
  ) +
  theme_bw() +
  theme(legend.position = "none") + 
  scale_color_manual(values = cols) +
  scale_y_log10() +
  scale_x_log10() +
  labs(
    x = "Tested | logarihmic scale", 
    y = "Positive | logarihmic scale",
    title = "Positive vs. Tested") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 



