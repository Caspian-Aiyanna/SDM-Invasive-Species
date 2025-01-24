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
#bioclimatic variables----
raster_files <- list.files("~/SDM/Arun/Arun/Data/worldclim/Current", 
                           pattern = "\\.tif$", 
                           full.names = TRUE)
envi_stack <- raster::stack(raster_files)
#ind<- st_read("~/SDM/Arun/Arun/Data/India Shape/india_st.shp")

pj <- read.csv("~/SDM/Arun/Arun/Data/Species_Occurrences/PJ.csv")
lc<-read.csv("~/SDM/Arun/Arun/Data/Species_Occurrences/LC.csv")
mc<-read.csv("~/SDM/Arun/Arun/Data/Species_Occurrences/MC.csv")


scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")


ESDM_Prosopis_juliflora <- ensemble_modelling(
  c('SVM', 'MARS', 'RF', 'GLM', 'ANN', 'CTA', 'GBM'),
  pj,
  envi_stack,
  Xcol = 'longitude',
  Ycol = 'latitude',
  pcol = NULL,
  rep = 10,
  cores = 10,
  cv = "holdout",
  cv.param = c(0.75, 1),
  ensemble.thresh = 0,
  verbose = FALSE
)
writeRaster(ESDM_Prosopis_juliflora@projection, filename = "~/SDM/Arun/Arun/Results/Prosopis_juliflora", format = "GTiff")
save(ESDM_Prosopis_juliflora, file = "Arun/Results/Prosopis_juliflora.rda")

ESDM_Muntingia_calabura <- ensemble_modelling(
  c('SVM', 'MARS', 'RF', 'GLM', 'ANN', 'CTA', 'GBM'),
  mc,
  envi_stack,
  Xcol = 'longitude',
  Ycol = 'latitude',
  pcol = NULL,
  rep = 10,
  cores = 10,
  cv = "holdout",
  cv.param = c(0.75, 1),
  ensemble.thresh = 0,
  verbose = FALSE
)
writeRaster(ESDM_Muntingia_calabura@projection, filename = "~/SDM/Arun/Arun/Results/Muntingia_calabura", format = "GTiff")
save(ESDM_Muntingia_calabura, file = "Arun/Results/Muntingia_calabura.rda")

ESDM_Lantana_camara <- ensemble_modelling(
  c('SVM', 'MARS', 'RF', 'GLM', 'ANN', 'CTA', 'GBM'),
  lc,
  envi_stack,
  Xcol = 'longitude',
  Ycol = 'latitude',
  pcol = NULL,
  rep = 10,
  cores = 10,
  cv = "holdout",
  cv.param = c(0.75, 1),
  ensemble.thresh = 0,
  verbose = FALSE
)
writeRaster(ESDM_Lantana_camara@projection, filename = "~/SDM/Arun/Arun/Results/Lantana_camara", format = "GTiff")
save(ESDM_Lantana_camara, file = "Arun/Results/Lantana_camara.rda")

load("Arun/Results/Lantana_camara.rda")
load("Arun/Results/Muntingia_calabura.rda")
load("Arun/Results/Prosopis_juliflora.rda")

# Create model_objects list
model_objects <- list(
  "Lantana_camara" = ESDM_Lantana_camara,
  "Muntingia_calabura" = ESDM_Muntingia_calabura,
  "Prosopis_juliflora" = ESDM_Prosopis_juliflora
)
# Load future climate data for each scenario
future_126 <- raster::stack(list.files("Arun/Data/ssp126", pattern = "\\.tif$", full.names = TRUE))
future_245 <- raster::stack(list.files("Arun/Data/ssp245", pattern = "\\.tif$", full.names = TRUE))
future_370 <- raster::stack(list.files("Arun/Data/ssp370", pattern = "\\.tif$", full.names = TRUE))
future_585 <- raster::stack(list.files("Arun/Data/ssp585", pattern = "\\.tif$", full.names = TRUE))

# Create future_stacks list
future_stacks <- list(
  future_126 = future_126,
  future_245 = future_245,
  future_370 = future_370,
  future_585 = future_585
)
names(future_stacks)

# Projection loop
for (species in names(model_objects)) {
  for (scen in scenarios) {
    scen_key <- paste0("future_", sub("ssp", "", scen))  # Remove 'ssp' from scenario name
    if (!scen_key %in% names(future_stacks)) {
      cat(sprintf("Environmental stack missing for %s - %s\n", species, scen))
      next  # Skip to the next iteration if the stack is missing
    }
    
    cat(sprintf("Starting projection for %s - %s\n", species, scen))
    
    future_proj <- SSDM::project(
      obj = model_objects[[species]],
      Env = future_stacks[[scen_key]],  # Use the corrected key
      uncertainty = TRUE,
      output.format = "model",
      SDM.projections = TRUE,
      cores =20,
      minimal.memory = FALSE,
      tmp = FALSE
    )
    
    # Save projection
    save(future_proj, 
         file = file.path("D:/Harin/Projects/SDM/Arun/Arun/Results", paste0(species, "_future_", scen, ".rda")))
    
    # Save raster
    writeRaster(future_proj@projection,
                filename = file.path("D:/Harin/Projects/SDM/Arun/Arun/Results/rasters", paste0(species, "_future_", scen)),
                format = "GTiff",
                overwrite = TRUE)
    
    cat(sprintf("Completed %s projection for %s\n", species, scen))
  }
}
