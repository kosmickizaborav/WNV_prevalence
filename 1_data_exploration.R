library(data.table)
library(ggplot2)
library(colorspace)
library(MetBrewer)
source("0_helper_functions.R")

data_dir <- here::here("Data")
graphs_dir <- file.path(data_dir, "Graphs")
dir.create(graphs_dir, showWarnings = F)

dt <- readxl::read_xlsx(
  file.path(data_dir, "WNV_Host_Susceptability_Review.xlsx"), 
  sheet = "observations") |> 
  janitor::clean_names()

dt <- as.data.table(dt)[class == "Aves"][bird_life_name != "NA"]
dt <- dt[ , .(
  bl_name = bird_life_name, total_tested, positive_m1 = positive_individuals_m1, 
  positive_m2 = positive_individuals_m2)]

# manual change for the species that cannot be found in the list of synonyms, 
# to incororate in the list later. 
# synonyms checked in avibase, but also found at birdlife datazone
#  Amazilia fimbriata = Chionomesa fimbriata (avibase, also at birdlife datazone)
#  Regulus calendula = Corthylio calendula (avibase, also at birdlife datazone)
#  Bonasa bonasia = Tetrastes bonasia
#  Amazilia tobaci = Saucerottia tobaci
#  Antigone canadensis = Grus canadensis
#  Antigone vipio = Grus vipio
dt[, bl_name := fcase(
  bl_name == "Amazilia fimbriata", "Chionomesa fimbriata",
  bl_name == "Regulus calendula", "Corthylio calendula",
  bl_name == "Bonasa bonasia", "Tetrastes bonasia",
  bl_name == "Amazilia tobaci", "Saucerottia tobaci",
  bl_name == "Antigone canadensis", "Grus canadensis",
  bl_name == "Antigone vipio", "Grus vipio",
  default = bl_name)]

# rename to new birdlife and add phylogeny
dt <- rename_to_birdlife(dt, species_name = "bl_name")[!is.na(birdlife_name)]
dt <- add_birdlife_phylogeny(dt, species_name = "birdlife_name")



# 1 - Add clade for color----------------------------------------------------

clades_all <- list(
  Paleognathae      = c("Struthioniformes","Rheiformes","Casuariiformes",
                        "Apterygiformes","Tinamiformes"),
  Galloanserae      = c("Galliformes","Anseriformes"),
  Strisores         = c("Caprimulgiformes","Podargiformes","Nyctibiiformes",
                        "Steatornithiformes","Apodiformes"),
  Columbaves        = c("Columbiformes","Pterocliformes","Mesitornithiformes",
                        "Cuculiformes","Musophagiformes","Otidiformes"),
  Gruimorphae       = c("Gruiformes","Opisthocomiformes","Eurypygiformes"),
  Aequorlitornithes = c("Charadriiformes","Phoenicopteriformes","Podicipediformes",
                        "Gaviiformes","Sphenisciformes","Procellariiformes",
                        "Phaethontiformes","Ciconiiformes","Suliformes",
                        "Pelecaniformes"),
  Telluraves        = c("Accipitriformes","Cathartiformes","Strigiformes",
                        "Coliiformes","Leptosomiformes", "Trogoniformes",
                        "Bucerotiformes","Coraciiformes","Piciformes",
                        "Passeriformes", "Psittaciformes", "Falconiformes")
)

# Assuming clades is your list as defined above
clades_all <- rbindlist(
  lapply(names(clades_all), function(clade) {
    data.table(order = clades_all[[clade]], clade = clade)
  }))

dt <- merge(dt, clades_all, by = "order", all.x = TRUE)
rm(clades_all)

clades_plot <- unique(dt[, .(order, clade)])

clade_count <- clades_plot[, .N, by = clade][order(-N)]

base_cols <- met.brewer("Archambault")
base_cols[3] <- "#6f9969"
base_cols[5] <- "#660d20"
names(base_cols) <- unique(clades_plot$clade)

cols_list <- list(
  Telluraves = met.brewer("Isfahan1", 11)[1:11], 
  Aequorlitornithes = met.brewer("Cassatt1", 10)[1:10], 
  Columbaves = met.brewer("VanGogh3", 5)[1:5],
  Galloanserae = met.brewer("Homer1")[c(4,5)],
  Gruimorphae = met.brewer("Homer1")[c(1,2)],
  Strisores = met.brewer("Kandinsky")[3],
  Paleognathae = met.brewer("Kandinsky")[4])

# Generate a named color vector for all orders
order_cols <- unlist(lapply(names(cols_list), function(cl) {
  orders <- clades_plot[clade == cl, order]
  cols <- cols_list[[cl]]
  names(cols) <- orders
  cols
}))


# plot 1 ------------------------------------------------------------------

# Summarize: number of species per family within each order
dt[, .(n_species = uniqueN(birdlife_name)), by = .(order, clade)] |> 
  ggplot(
    aes(x = reorder(order, -n_species), 
        y = n_species, group = clade, fill = order)
  ) +
  geom_bar(stat = "identity", color = "gray33") +
  geom_text(aes(label = n_species), vjust = -0.5) +
  labs(
    x = "Order",
    y = "Number of Species",
    title = "1426 bird species from 129 families and 32 orders",
  ) +
  scale_fill_manual(values = order_cols) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold"), 
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) 

ggsave(
  filename = file.path(graphs_dir, "1_species_counts_by_order.png"),
  width = 10, height = 6, dpi = 300, bg = "transparent"
)


library(patchwork)  # For arranging multiple plots

# Assume your data.table is called dt, with columns birdlife_name (species), family, order

# Count unique species per family within each order
sp_per_family <- dt[
  , .(n_species = uniqueN(birdlife_name)), by = .(order, family)][
  , ':=' (
    total_species = sum(n_species),
    n_family = uniqueN(family)), by = order][total_species > 10][
  , order_family := paste(order, family, sep = "_")][
  , order_lab := paste0(order, " [N = ", n_family, "]")]

get_fam_pal <- function(base_color, n) {
  lighten(base_color, seq(0, 0.4, length.out = n))
}

# Create a color vector: names are "order_family", values are hex colors
family_cols_list <- sp_per_family[
  , {
    fams <- unique(family)
    cols <- get_fam_pal(order_cols[unique(order)], length(fams))
    setNames(cols, paste(unique(order), fams, sep = "_"))
  },
  by = order
]
# If the above gives you a list, flatten it:
family_cols <- unlist(family_cols_list$V1)

ggplot(sp_per_family, aes(x = "", y = n_species, fill = order_family)) +
  geom_bar(stat = "identity", width = 1, color = "gray33") +
  coord_polar("y") +
  facet_wrap(~order_lab, scales = "free") +
  scale_fill_manual(values = family_cols) + 
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "", title = "Number of species per family and order") 

ggsave(
  filename = file.path(graphs_dir, "1_species_per_family_pie.png"),
  width = 10, height = 10, dpi = 300, bg = "transparent"
)


cols_dt <- rbindlist(lapply(unique(sp_per_family$order), function(ord) {
  fams <- unique(sp_per_family[order == ord, family])
  cols <- get_fam_pal(order_cols[ord], length(fams))
  data.table(
    order = ord, 
    family = fams, 
    color = cols,
    order_family = paste(ord, fams, sep = "_")
  )}), fill = T)

# numbers -----------------------------------------------------------------

cols <- cols_dt[, color]
names(cols) <- cols_dt$family

copy(dt)[, positive_m1 := ifelse(positive_m1 == 0, 0.1, positive_m1)] 

copy(dt)[, positive_m1 := ifelse(positive_m1 == 0, 0.1, positive_m1)] |> 
  ggplot() +
  geom_point(
    aes(x = total_tested, y = positive_m1, color = as.factor(family)), alpha = 0.5
  ) +
  theme_bw() +
  theme(legend.position = "none") + 
  scale_color_manual(values = cols) + 
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



# traits ------------------------------------------------------------------

traits <- fread(file.path(data_dir, "Traits_filtered.csv")) |> 
  janitor::clean_names()

dt |> 
  ggplot() +
  geom_point(
    aes(x = total_tested, y = positive_m1, color = as.factor(family)), alpha = 0.5
  ) +
  theme_bw() +
  theme(legend.position = "none") + 
  scale_color_manual(values = cols) +
  coord_cartesian(xlim = c(0, 1000)) 

"#4DB5C5"


prev <- readxl::read_xlsx("Data/Palearctic_prevalence bird_species.xlsx") |> 
  janitor::clean_names()

prev <- as.data.table(prev)
prev <- prev[, .(bl_name = bird_life_name, prevalence = group_prevalence)] 

# prev[, bl_name := fcase(
#   bl_name == "Amazilia fimbriata", "Chionomesa fimbriata",
#   bl_name == "Regulus calendula", "Corthylio calendula",
#   bl_name == "Bonasa bonasia", "Tetrastes bonasia",
#   bl_name == "Amazilia tobaci", "Saucerottia tobaci",
#   bl_name == "Antigone canadensis", "Grus canadensis",
#   bl_name == "Antigone vipio", "Grus vipio",
#   default = bl_name)]


prev_dt <- dt[, .(
  positive = sum(positive_m1, na.rm = TRUE), 
  total = sum(total_tested, na.rm = TRUE), 
  prev = sum(positive_m1, na.rm = TRUE) / sum(total_tested, na.rm = TRUE)), 
              by = .(bl_name, family)][total >= 10]
prev_traits <- merge(prev_dt, traits, by.x = "bl_name", by.y = "scientific_name")




prev_traits |> 
  ggplot() + 
  geom_boxplot(aes(y = prev, x = habitat, fill = habitat), alpha = 0.8) +
  theme_bw() +
  theme(legend.position = "none") +
  MetBrewer::scale_fill_met_d("Renoir") +
  labs(x = "habitat", y = "positive/tested summarized by species")

ggsave(
  filename = file.path(graphs_dir, "1_prev_habitat_boxplot.png"),
  width = 11, height = 6, dpi = 300, bg = "transparent"
)

prev_traits |> 
  ggplot() + 
  geom_boxplot(aes(y = prev, x = primary_lifestyle, fill = primary_lifestyle), alpha = 0.8) +
  theme_bw() +
  theme(legend.position = "none") +
  MetBrewer::scale_fill_met_d("Renoir") +
  labs(x = "primary lifestyle", y = "positive/tested summarized by species")

ggsave(
  filename = file.path(graphs_dir, "1_prev_primlife_boxplot.png"),
  width = 8, height = 6, dpi = 300, bg = "transparent"
)


prev_traits |> 
  ggplot() + 
  geom_boxplot(aes(y = prev, x = as.factor(migration), fill = as.factor(migration)), alpha = 0.8) +
  theme_bw() +
  theme(legend.position = "none") +
  MetBrewer::scale_fill_met_d("Renoir") +
  labs(x = "territoriality", y = "positive/tested summarized by species")

ggsave(
  filename = file.path(graphs_dir, "1_prev_forage_boxplot.png"),
  width = 13, height = 6, dpi = 300, bg = "transparent"
)


prev_traits |> 
  ggplot() + 
  geom_point(aes(y = prev, x = abundance_estimate, color = family), alpha = 0.8) +
  theme_bw() +
  scale_color_manual(values = cols) +
  scale_x_log10() +
  theme(legend.position = "none") +
  labs(x = "maximum longevity", y = "positive/tested summarized by species")

ggsave(
  filename = file.path(graphs_dir, "1_prev_long_points.png"),
  width = 6, height = 6, dpi = 300, bg = "transparent"
)

dt <- rename_to_birdlife(dt, species_name = "bl_name")[!is.na(birdlife_name)]
