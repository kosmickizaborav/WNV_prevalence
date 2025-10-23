#' **FUNCTIONS taken from the project WNV_dispersal!**


# Preparing birdlife classification ---------------------------------------

# There were some species missing, so added them manually, used the code that
# was taken from WNV_dispersal project
# synonyms added = synonym, birdlife_name: 
# "Amazilia fimbriata", "Chionomesa fimbriata",
# "Regulus calendula", "Corthylio calendula",
# "Bonasa bonasia", "Tetrastes bonasia",
# "Amazilia tobaci", "Saucerottia tobaci",
# "Antigone canadensis", "Grus canadensis",
# "Antigone vipio", "Grus vipio",

if(!file.exists(here::here("Data", "00_birdlife_classification.csv"))){
  
  library(tidyverse)
  
  # BIRDLIFE classification downloaded from: 
  # https://datazone.birdlife.org/about-our-science/taxonomy
  
  bird_file <- "Handbook of the Birds of the World and BirdLife International Digital Checklist of the Birds of the World_Version_9.xlsx"
  
  birdlife <- file.path(here::here("Data", bird_file)) |>
    readxl::read_xlsx(skip = 2) |>
    # take only the species level, not subspecies
    janitor::clean_names() |>
    rename(sp_status = x2024_iucn_red_list_category) |>
    select(seq, scientific_name, synonyms, family_name, order, sp_status) |>
    filter(!is.na(seq)) |>
    distinct(seq, .keep_all = T) |>
    rename(family = family_name, birdlife_seq = seq) |>
    # remove information in parenthesis is added to the synonyms column
    mutate(synonyms = str_remove_all(synonyms, "\\([^\\)]*\\)")) |>
    # change all the punctuation to ";"
    mutate(synonyms = str_replace_all(synonyms, "[:punct:]", ";")) |> 
    # if the species had synonym already, added them to the list
    mutate(
      synonyms = case_when(
        scientific_name == "Curruca curruca" ~ str_c(
          synonyms, ";Sylvia curruca", sep = ";"),
        scientific_name == "Chionomesa fimbriata" ~ str_c(
          synonyms, "Amazilia fimbriata", sep = ";"),
        .default = synonyms
      )
    ) |>
    # extract synonyms one per row
    separate_longer_delim(synonyms, delim = ";") |>
    # remove labels for variants and subspecies
    # remove synonyms that are not a full species name but just an epitaph or
    # genus, it's the product of previous steps
    mutate(synonyms = str_remove_all(synonyms, "\\s(sp$|var$|ssp$)")) |>
    # making sure there is no extra white spaces
    mutate(across(everything(), str_squish)) |>
    mutate(
      # removing the synonyms that are just the epitaph or genus name
      # or the ones that are identical to the scientific name
      synonyms = case_when(
        str_count(synonyms, " ") !=1 ~ NA,
        scientific_name == synonyms ~ NA,
        .default = synonyms
      ),
      # change the status to NR - not recognized and R - regonized
      sp_status = if_else(sp_status == "NR", "NR", "R"),
      order = str_to_title(order)
    ) |>
    distinct(scientific_name, synonyms, family, sp_status, .keep_all = T) |> 
    mutate(
      # checked that they don't have any synonyms, so I just manually added them
      synonyms = case_when(
        scientific_name == "Myiopsitta monachus" ~ "Psittacula monachus",
        scientific_name == "Cyanistes caeruleus" ~ "Parus caeruleus",
        scientific_name == "Larus melanocephalus" ~ "Ichthyaetus melanocephalus",
        scientific_name == "Larus genei" ~ "Chroicocephalus genei",
        scientific_name == "Larus audouinii" ~ "Ichthyaetus audouinii",
        scientific_name == "Corvus monedula" ~ "Coloeus monedula",
        scientific_name == "Corthylio calendula" ~ "Regulus calendula",
        scientific_name == "Tetrastes bonasia" ~ "Bonasa bonasia",
        scientific_name == "Saucerottia tobaci" ~ "Amazilia tobaci",
        scientific_name == "Grus canadensis" ~ "Antigone canadensis",
        scientific_name == "Grus vipio" ~ "Antigone vipio",
        .default = synonyms
      ),
      synonyms = if_else(
        str_detect(scientific_name, "Curruca") & is.na(synonyms),
        paste("Sylvia", str_split_i(scientific_name, " ", 2)),
        synonyms
      )
    ) |>
    # check how many times the species appears,
    # whether it is both listed as recognized and not
    # how many distinct synonyms it has
    mutate(
      n_sci = n(),
      both_r_nr = sum(c("R", "NR") %in% sp_status) == 2,
      dist_syn = n_distinct(synonyms, na.rm = T),
      .by = scientific_name
    ) |>
    # check unique combos
    mutate(n_pair = n(), .by = c(scientific_name, synonyms)) |>
    mutate(
      sci_in_syn = scientific_name %in% unique(synonyms),
      syn_in_sci = synonyms %in% unique(scientific_name)
    ) |>
    mutate(
      sci_name = case_when(
        # if scientific name is both listed as recognized and not recognized,
        # remove the non recognized if it has no additional synonyms
        n_sci > 1 & both_r_nr & sp_status == "NR" & is.na(synonyms) ~ NA,
        # if it does have synonyms, make sure they are identical as in recognized
        n_sci > 1 & both_r_nr & sp_status == "NR" & n_pair > 1 ~ NA,
        # if the scientific name is duplicated but one of the entries doesn't
        # have synonym provided, remove it
        n_sci > 1 & dist_syn > 0 & is.na(synonyms) ~ NA,
        # if the scientific name is provided also found in the synonyms, but
        # it's listed as not recognized, remove it if there is no synonym provided
        sci_in_syn & sp_status == "NR" & is.na(synonyms) ~ NA,
        .default = scientific_name
      )
    ) |>
    filter(!is.na(sci_name)) |>
    mutate(
      sp_status = if_else(sp_status == "NR" & both_r_nr == T, "R", sp_status),
      syn_in_sci = synonyms %in% unique(scientific_name),
      synonyms = if_else(syn_in_sci, NA, synonyms)
    ) |>
    filter(
      !(n() > 1 & n_distinct(synonyms, na.rm = T) > 0 & is.na(synonyms)),
      .by = scientific_name
    ) |>
    select(birdlife_seq, scientific_name, synonyms, family, order, sp_status) |>
    # assign a id number just for easier organization
    rename(synonym = synonyms)
  
  birdlife |>
    write_csv(here::here("Data", "00_birdlife_classification.csv"))

}


# FUNCTION: add_birdlife_phylogeny ---------------------------------------------

add_birdlife_phylogeny <- function(
    df = NULL, 
    species_name, 
    birdlife_file = here::here("Data", "00_birdlife_classification.csv")
) {
  
  # Validate inputs
  if (is.null(df) && !is.vector(species_name)) {
    stop("Either `df` must be provided or `species_name` must be a vector.")
  }
  
  if(!is.null(df) && !is.data.table(df)){
    stop("`df` must be a data.table or NULL.")
  }
  
  if (!is.null(df) && !species_name %in% names(df)) {
    stop("The specified `species_name` column does not exist in `df`.")
  }
  
  birdlife <- data.table::fread(birdlife_file)
  birdlife <- unique(birdlife[, .(scientific_name, family, order)])
    
  # If `df` is NULL and `species_name` is a vector, create a data.table
  if (is.null(df) && is.vector(species_name)) {
    df <- data.table(birdlife_name = species_name)
    species_name <- "birdlife_name"
  }
  
  
  # Merge with birdlife synonyms
  df <- merge(
    df,
    birdlife,
    by.x = species_name,
    by.y = "scientific_name",
    all.x = TRUE
  )
  
  return(df)
}


# FUNCTION: add_birdlife_phylo --------------------------------------------

rename_to_birdlife <- function(
    df = NULL, 
    species_name, 
    max_dist = 2, 
    birdlife_file = here::here("Data", "00_birdlife_classification.csv")
) {
  # Validate inputs
  if (is.null(df) && !is.vector(species_name)) {
    stop("Either `df` must be provided or `species_name` must be a vector.")
  }
  
  if(!is.null(df) && !is.data.table(df)){
    stop("`df` must be a data.table or NULL.")
  }
  
  if (!is.null(df) && !species_name %in% names(df)) {
    stop("The specified `species_name` column does not exist in `df`.")
  }
  
  birdlife <- fread(birdlife_file)
  
  # Prepare birdlife data in long format for matching
  bln <- melt(
    birdlife[, .(scientific_name, synonym, sp_status)],
    id.vars = c("sp_status"),
    measure.vars = c("scientific_name", "synonym"),
    variable.name = "name_type",
    value.name = "sci_name",
    na.rm = TRUE
  )[
    , name_type := ifelse(
      sp_status == "R",
      paste("recognized", gsub("_", " ", name_type)),
      paste("not recognized", gsub("_", " ", name_type))
    )
  ][!duplicated(sci_name)]
  
  # Matching function for approximate string matching
  match_names <- function(sp) {
    loc <- stringdist::amatch(
      sp, bln[, sci_name], maxDist = max_dist, matchNA = FALSE, method = "lv"
    )
    bln[loc, sci_name]
  }
  
  # If `df` is NULL and `species_name` is a vector, create a data.table
  if (is.null(df) && is.vector(species_name)) {
    df <- data.table(original_name = species_name)
    species_name <- "original_name"
  }
  
  # Add birdlife_name column and match with birdlife synonyms
  df[, birdlife_name := get(species_name)]
  df[, birdlife_name := fifelse(
    !birdlife_name %in% bln[, sci_name],
    match_names(birdlife_name),
    birdlife_name
  )]
  
  # Merge with birdlife synonyms
  df <- merge(
    df,
    bln[, .(sci_name, name_type)],
    by.x = "birdlife_name",
    by.y = "sci_name",
    all.x = TRUE
  )
  
  # Replace synonyms with official scientific names
  df[, birdlife_name := fifelse(
    grepl("synonym", name_type),
    birdlife$scientific_name[match(birdlife_name, birdlife$synonym)],
    birdlife_name
  )]
  
  return(df)
}



# match_to_avilist --------------------------------------------------------

match_to_avilist <- function(
    df = NULL, 
    species_name, 
    avilist_file = here::here("Data", "Avilist_match.xlsx"),
    phylo_file = here::here("Data", "AviList-v2025-11Jun-short.xlsx"),
    add_phylo = F
){
  
  # Validate inputs
  if (is.null(df) && !is.vector(species_name)) {
    stop("Either `df` must be provided or `species_name` must be a vector.")
  }
  
  if(!is.null(df) && !is.data.table(df)){
    stop("`df` must be a data.table or NULL.")
  }
  
  if (!is.null(df) && !species_name %in% names(df)) {
    stop("The specified `species_name` column does not exist in `df`.")
  }
  
  # If `df` is NULL and `species_name` is a vector, create a data.table
  if (is.null(df) && is.vector(species_name)) {
    df <- data.table(original_name = species_name)
    species_name <- "original_name"
  }
  
  avilist <- setDT(readxl::read_xlsx(avilist_file)) |> 
    janitor::clean_names()
  
  # set match_spp to NA if exactly the same as avilist_species
  avilist[species_avilist == match_spp, match_spp := NA_character_]
  
  # there are sometimes multiple identical match_spp leading to different 
  # avilist_species, to avoid that, we select the species that has the same 
  # epithet as in match_spp, if none match, just take the first one
  avilist[, avilist_epithet := sub("^[^ ]+ +", "", species_avilist)]
  avilist[, match_epithet := sub("^[^ ]+ +", "", match_spp)]
  
  # For duplicated match_spp (ignoring NA), keep only one match per group:
  avilist[!is.na(match_spp), 
     flag := {
       # If any epithet matches, keep it; else keep first only
       match_rows <- which(avilist_epithet == match_epithet)
       if (length(match_rows) > 0) {
         seq_len(.N) %in% match_rows
       } else {
         seq_len(.N) == 1
       }
     }, 
     by = match_spp]
  

  # set match_spp to NA where flag is FALSE
  avilist[!is.na(match_spp) & !flag, match_spp := NA_character_]
  
  # remove helper columns
  avilist[, c("avilist_epithet", "match_epithet", "flag") := NULL]
  
  df[, match_spp := get(species_name)]
  df[, name_type := fcase(
    match_spp %in% avilist$species_avilist, "avilist", 
    match_spp %in% avilist$match_spp, "match", 
    default = NA_character_)]
  
  df <- merge(
    df, 
    avilist[, .(species_avilist, match_spp)],
    by = "match_spp", 
    all.x = T
  )
  
  df[, avilist_name := fcase(
    name_type == "avilist", match_spp, 
    name_type == "match", species_avilist, 
    default = NA_character_)]
  
  df[, c("species_avilist", "match_spp") := NULL]
  
  if(add_phylo == T){
    
    phylo_dt <- readxl::read_xlsx(phylo_file) |> 
      janitor::clean_names()
    
    phylo_dt <- setDT(phylo_dt)[
      , .(avilist_name = scientific_name, family, order)]
    
    df <- merge(df, phylo_dt, by = "avilist_name", all.x = T)
    
  }
  
  return(df)
}



