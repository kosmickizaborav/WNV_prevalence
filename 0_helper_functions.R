#' **FUNCTIONS taken from the project WNV_dispersal!**

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
    max_dist = 1, 
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


