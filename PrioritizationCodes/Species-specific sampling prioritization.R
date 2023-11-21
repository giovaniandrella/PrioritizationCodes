##%######################################################%##
#                                                          #
####     Species-specific sampling prioritization       ####
#                                                          #
##%######################################################%##

{
  require(terra)
  require(dplyr)
  require(spatstat)
  require(data.table)
  require(spatstat.explore)
  require(spatstat.geom)
}

# Set working directory to PrioritizationCodes folder
setwd("C:/Users/my_computer/.../.../PrioritizationCodes")

# function to standardize rasters between 0 and 1
stand_r <- function(r) {
  minn <- terra::global(r, min, na.rm = TRUE)[1, 1]
  maxx <- terra::global(r, max, na.rm = TRUE)[1, 1]
  r <- ((r - minn) / (maxx - minn))
  return(r)
}



##%######################################################%##
#                                                          #
####     Calculate point density for each species       ####
#                                                          #
##%######################################################%##
# Process all species based on occurrence density

# Read occurrences
occ <- data.table::fread("./occurrence_data.txt") %>%
  dplyr::tibble()
occ

# Read basic layer
basic_l <- "./basic_l.tif" %>% terra::rast()
plot(basic_l)

# List of species names
spp <- occ$Species %>%
  unique() %>%
  sort()

# Create directory for saving raster with species destiny
dir.create("occ_dens_sp")

# Loop for processing each species

for (i in 1:length(spp)) {
  message("Processing species ", i)
  occ_spp <- occ %>%
    dplyr::filter(Species == spp[i]) %>%
    unique()
  
  # Create point pattern object
  ppp_poinst <-
    spatstat.geom::ppp(
      x = occ_spp$Longitude,
      y = occ_spp$Latitude,
      window = as.owin(c(ext(basic_l)[1:2], ext(basic_l)[3:4]))
    )

  # Cross validated bandwidth selection for kernel density
  if (nrow(occ_spp) > 1) {
    hmax <- spatstat.explore::rmax.rule("K", W = spatstat.geom::Window(ppp_poinst),
      lambda = spatstat.geom::npoints(ppp_poinst) / spatstat.geom::area(Window(ppp_poinst))
    )
    sgm <- spatstat.explore::bw.diggle(X = ppp_poinst, hmax = hmax)
  } else {
    sgm <- 0.02
  }

  # Calculate density
  p_density <-
    spatstat.explore::density.ppp(
      ppp_poinst,
      sigma = sgm,
      edge = TRUE,
      at = "pixel"
    ) # check how set output resolution
  p_density <- terra::rast(p_density)
  p_density <- p_density %>% stand_r()
  terra::crs(p_density) <- terra::crs(basic_l)
  p_density <- p_density %>%
    terra::resample(., basic_l) %>%
    terra::mask(., basic_l)
  p_density <- 1 - p_density # invert density
  plot(p_density)
  writeRaster(p_density, file.path("./occ_dens_sp", paste0(spp[i], ".tif")), overwrite = TRUE)
}


## %######################################################%##
#                                                          #
####  Prioritization (suitability * density of points)  ####
#                                                          #
## %######################################################%##

# Create directory to save prioritization for each species
dir.create("prioritization_sp")

# species suitability
suit <- "./models/" %>% 
  list.files(., full.names = TRUE, pattern = ".tif$")
names(suit) <- basename(suit) %>% 
  gsub(".tif$", "", .) %>% 
  gsub("_", " ", .)

# species density
dens <- "./occ_dens_sp/" %>% 
  list.files(., full.names = TRUE, pattern = ".tif$")
names(dens) <- basename(dens) %>% 
  gsub(".tif$", "", .) %>% 
  gsub("_", " ", .)

length(suit) # Number of distribution rasters
length(dens) # Number of point density rasters

# Calculate prioritizatio for each species 
for (i in 1:length(suit)) {
  message("Processing species ", i)
  s <- terra::rast(suit[i])
  d <- terra::rast(dens[i])
  sdr <- s * d
  writeRaster(sdr, file.path("prioritization_sp", paste0(names(suit[i]), ".tif")), overwrite = TRUE)
}


# Sum of species prioritized maps 
priort <- "./prioritization_sp" %>% 
  list.files(., full.names = TRUE, pattern = ".tif$")
names(priort) <- basename(priort) %>% 
  gsub(".tif$", "", .) %>% 
  gsub("_", " ", .)

prio_all_sp <- terra::rast(priort[1])
for (i in 2:length(priort)) {
  prio_all_sp <- prio_all_sp + terra::rast(priort[i])
}

prio_all_sp <- stand_r(prio_all_sp)
plot(prio_all_sp, main="Prioritization map")

# Create directory to save final results
dir.create("results")
writeRaster(prio_all_sp, "./results/Sample_priority_map.tif")



##%######################################################%##
#                                                          #
####     Adding constraints to prioritization maps      ####
#                                                          #
##%######################################################%##

rdist <- "./constraints/distance_roads.tif" %>% terra::rast()
lcover <- "./constraints/land_cover.tif" %>% terra::rast()
pas <- "./constraints/protected_areas.tif" %>% terra::rast()

# Priority sampling map
prio_all_sp <- "./results/Sample_priority_map.tif" %>% terra::rast()

# Distance to road
plot(prio_all_sp, main = "Unconstrained prioritization map")
plot(rdist, main = "Distance to roads")
rdist_bin <- rdist < 5
plot(rdist_bin, main = "Distance to roads < 5km")
plot(prio_all_sp * rdist_bin, main = "Prioritization map constrained by roads")

# Land cover
plot(lcover, main ="Land cover") # Higher values represent more conserved landscapes 
lcover_bin <- lcover > 0.8
plot(lcover_bin,  main ="Land cover >80% of lanscape conserved")
plot(prio_all_sp * lcover_bin, main = "Prioritization map constrained by land use") 

# Protected Areas
plot(pas, main= "Protected Areas") # Lower values represent no unprotected cells
pas_inv <- 1 - pas # invert values
plot(pas_inv, main= "Unprotected Areas") # Lower values represent no unprotected cells
plot(prio_all_sp * pas_inv, main = "Prioritization map constrained by unprotected areas")

# Prioritization map combining all constraints
prio_const <- (prio_all_sp * rdist_bin * lcover_bin * pas_inv)
plot(prio_const)
terra::writeRaster(prio_const, "./results/constrained_sampling_priority_map.tif")
