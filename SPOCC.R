library(data.table)
library(dplyr)
library(terra)
library(spocc)
library(sf)
library(ggplot2)
library(readr)
library(rnaturalearth)
library(ggplot2)
library(SSDM)
library(dismo)

setwd("~/SDM/Arun")
#species occurence----
species_list <- c("Lantana camara", "Prosopis juliflora", "Muntingia calabura")
output_dir <- "Arun/Data/Species_Occurrences"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ind<- st_read("~/SDM/Arun/Arun/Data/India Shape/india_st.shp")
india_fixed <- st_make_valid(ind)
india_buffer <- st_buffer(india_fixed, 0)
india_geom <- st_as_text(st_convex_hull(st_union(india_buffer)))

#Function to fetch species data
fetch_species_data <- function(species_name, geometry) {
  cat("Fetching data for", species_name, "\n")
  
  occ_data <- spocc::occ(query = species_name, 
                         geometry = geometry, 
                         from = c('gbif', 'inat', 'ebird', 'OBIS', 'iDigBio'), 
                         limit = 2500000)
  
  occ_df <- spocc::occ2df(occ_data)
  
  if (nrow(occ_df) == 0) {
    return(data.frame(Species = character(),
                      longitude = numeric(),
                      latitude = numeric()))
  }
  
  # Convert to data.table
  filtered_dt <- as.data.table(occ_df)
  filtered_dt[, year := as.numeric(substr(as.character(date), 1, 4))]
  filtered_dt[, c("longitude", "latitude") := .(as.numeric(longitude), as.numeric(latitude))]
  filtered_dt <- filtered_dt[year >= 1975 & year <= 2024]
  filtered_dt <- filtered_dt[, .(Species = name, longitude, latitude)]
  
  return(filtered_dt)
}

#Process each species
for (species in species_list) {
  species_data <- fetch_species_data(species, india_geom)
  output_file <- file.path(output_dir, paste0(gsub(" ", "_", species), "_occurences.csv"))
  write_csv(species_data, output_file)
  cat("Data for", species, "saved to", output_file, "\n")
}

#bioclimatic variables----
raster_files <- list.files("~/SDM/Arun/Arun/Data/worldclim/Current", 
                           pattern = "\\.tif$", 
                           full.names = TRUE)
#1st model for Prosopis juliflora----
species_data <- read.csv("~/SDM/Arun/Arun/Data/Species_Occurrences/Prosopis_juliflora_occurences.csv")
envi_stack <- raster::stack(raster_files)

ESDM_Prosopis_juliflora <- ensemble_modelling(
  c('RF'),
  species_data,
  envi_stack,
  Xcol = 'Longitude',
  Ycol = 'Latitude',
  pcol = NULL,
  name = Species,
  rep = 10,
  cores = 10,
  cv = "holdout",
  cv.param = c(0.75, 1),
  ensemble.thresh = 0,
  verbose = FALSE
)
writeRaster(ESDM_Prosopis_juliflora@projection, filename = "~/SDM/Arun/Arun/Results/Prosopis_juliflora", format = "GTiff")
save(ESDM_Prosopis_juliflora, file = "Arun/Results/Prosopis_juliflora.rda")

#model building in loop----
  # Create a list of your species
  # Setup paths and list
species_list <- c("Lantana_camara", "Muntingia_calabura")
data_path <- "Arun/Data/Species_Occurrences/"
results_path <- "Results/"
raster_path <- "Results/rasters/"

# Create directories if they don't exist
dir.create(results_path, showWarnings = FALSE)
dir.create(raster_path, showWarnings = FALSE)

# Process each species
for (species in species_list) {
  
  # Run ensemble modeling
  ESDM_present <- ensemble_modelling(
    c('RF'),
    model_data,
    envi_stack,
    Xcol = 'Longitude',
    Ycol = 'Latitude',
    pcol = NULL,
    name = Species,
    rep = 10,
    cores = 10,
    cv = "holdout",
    cv.param = c(0.75, 1),
    ensemble.thresh = 0,
    verbose = FALSE
  )
  
  # Save results
  save(ESDM_present, 
       file = file.path(results_path, paste0(species, "_present.rda")))
  
  writeRaster(ESDM_present@projection, 
              filename = file.path(raster_path, paste0(species, "_present")), 
              format = "GTiff",
              overwrite = TRUE)
  
  # Print progress
  cat(sprintf("Completed processing for %s\n", species))
}



#future model----
#for Prosopis juliflora----
base_path <- "~/SDM/Arun/Arun"
scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")
# Load future climate data for each scenario
future_126 <- raster::stack(list.files("~/SDM/Arun/Arun/Data/ssp126", pattern = "\\.tif$", full.names = TRUE))
future_245 <- raster::stack(list.files("~/SDM/Arun/Arun/Data/ssp245", pattern = "\\.tif$", full.names = TRUE))
future_370 <- raster::stack(list.files("~/SDM/Arun/Arun/Data/ssp370", pattern = "\\.tif$", full.names = TRUE))
future_585 <- raster::stack(list.files("~/SDM/Arun/Arun/Data/ssp585", pattern = "\\.tif$", full.names = TRUE))
# # Rename layers to match the training data names
# names(future_126) <- c("india_bioclim_1", "india_bioclim_2", "india_bioclim_3", "india_bioclim_4", 
#                        "india_bioclim_5", "india_bioclim_6", "india_bioclim_7", "india_bioclim_8", 
#                        "india_bioclim_9", "india_bioclim_10", "india_bioclim_11", "india_bioclim_12", 
#                        "india_bioclim_13", "india_bioclim_14", "india_bioclim_15", "india_bioclim_16", 
#                        "india_bioclim_17", "india_bioclim_18", "india_bioclim_19")
# 
# names(future_245) <- names(future_126)
# names(future_370) <- names(future_126)
# names(future_585) <- names(future_126)
# 
# future_stacks <- list(
#   "future_126" = future_126,
#   "future_245" = future_245,
#   "future_370" = future_370,
#   "future_585" = future_585
# )

scenarios <- c("126", "245", "370", "585")

for (scen in scenarios) {
  future_proj <- SSDM::project(
    obj = ESDM_Prosopis_juliflora,
    Env = future_stacks[[paste0("future_", scen)]],
    uncertainty = TRUE,
    output.format = "model",
    SDM.projections = TRUE,
    cores = 16,  # Adjusted to match your CPU capacity
    minimal.memory = FALSE,
    tmp = FALSE
  )
  
  save(future_proj, 
       file = file.path("Results", 
                        paste0("Prosopis_juliflora_future_", scen, ".rda")))
  
  writeRaster(future_proj@projection,
              filename = file.path("Results", 
                                   paste0("Prosopis_juliflora_future_", scen)),
              format = "GTiff",
              overwrite = TRUE)
  
  cat(sprintf("Completed Prosopis_juliflora projection for ssp%s\n", scen))
}



#for rest of the species----

# Load all models
load("~/SDM/Arun/Arun/Results/Lantana_camara_present.rda")
load("~/SDM/Arun/Arun/Results/Muntingia_calabura_present.rda")

# Define species list and scenarios
species_list <- c("Lantana_camara", "Muntingia_calabura")
scenarios <- c("126", "245", "370", "585")

#named list of the loaded models
model_objects <- list(
  "Lantana_camara" = ESDM_present,
  "Muntingia_calabura" = ESDM_present,
)

#projection loop
for (species in names(model_objects)) {
  for (scen in scenarios) {
    future_proj <- SSDM::project(
      obj = model_objects[[species]],
      Env = future_stacks[[paste0("future_", scen)]],
      uncertainty = TRUE,
      output.format = "model",
      SDM.projections = TRUE,
      cores = 10,
      minimal.memory = FALSE,
      tmp = FALSE
    )
    
    # Saveprojection
    save(future_proj, 
         file = file.path("Results", paste0(species, "_future_ssp", scen, ".rda")))
    
    # Saveraster
    writeRaster(future_proj@projection,
                filename = file.path("Results", paste0(species, "_future_ssp", scen)),
                format = "GTiff",
                overwrite = TRUE)
    
    cat(sprintf("Completed %s projection for ssp%s\n", species, scen))
  }
}

#STACKED SPECIES DISTRIBUTION MODEL
