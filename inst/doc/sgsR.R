## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----echo=FALSE---------------------------------------------------------------

p <- c("`mraster`", "`sraster`", "`access`", "`existing`", "plot")

d <- c("Metric raster(s)", "Stratified raster", "Linear vectors representing access routes", "Existing sample units", "Visually displays raster and samples")

df <- data.frame(Parameter = p, Description = d)

knitr::kable(df, align = 'c')


## ----warning=F,message=F------------------------------------------------------
library(sgsR)
library(terra)
library(sf)

#--- Load mraster from internal data ---#
r <- system.file("extdata", "mraster.tif", package = "sgsR")

#--- load mraster using the terra package ---#
mraster <- terra::rast(r)

## ----warning=F,message=F------------------------------------------------------
#--- apply kmeans algorithm to metrics raster ---#
sraster <- strat_quantiles(mraster = mraster$zq90, # use mraster as input for sampling
                           nStrata = 4, # algorithm will produce 4 strata
                           plot = TRUE) # algorithm will plot output


## ----warning=F,message=F------------------------------------------------------
#--- apply stratified sampling ---#
existing <- sample_strat(sraster = sraster, # use mraster as input for sampling
                         nSamp = 200, # request 200 samples be taken
                         mindist = 100, # define that samples must be 100 m apart
                         plot = TRUE) # algorithm will plot output


## ----warning=F,message=F------------------------------------------------------
a <- system.file("extdata", "access.shp", package = "sgsR")

#--- load the access vector using the sf package ---#
access <- sf::st_read(a)

## ----warning=F,message=F------------------------------------------------------
terra::plot(mraster$zq90)
terra::plot(access, add = TRUE, col = "black")

## ----pipe, eval= FALSE--------------------------------------------------------
#  #--- non piped ---#
#  sraster <- strat_quantiles(mraster = mraster$zq90, # use mraster as input for sampling
#                             nStrata = 4) # algorithm will produce 4 strata
#  
#  existing <- sample_strat(sraster = sraster, # use mraster as input for sampling
#                           nSamp = 200, # request 200 samples be taken
#                           mindist = 100) # define that samples must be 100 m apart
#  
#  extract_metrics(mraster = mraster,
#                  existing = existing)
#  
#  #--- piped ---#
#  strat_quantiles(mraster = mraster$zq90, nStrata = 4) %>%
#    sample_strat(., nSamp = 200, mindist = 100) %>%
#    extract_metrics(mraster = mraster, existing = .)
#  

