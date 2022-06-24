## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=F,message=F,echo=FALSE,results=FALSE-----------------------------
library(sgsR)
library(terra)

#--- Load mraster and access files ---#
r <- system.file("extdata", "mraster.tif", package = "sgsR")

#--- load the mraster using the terra package ---#
mraster <- terra::rast(r)

a <- system.file("extdata", "access.shp", package = "sgsR")

#--- load the access vector using the sf package ---#
access <- sf::st_read(a, quiet = TRUE)

#--- apply quantiles algorithm to metrics raster ---#
sraster <- strat_quantiles(mraster = mraster$zq90, # use mraster as input for sampling
                           nStrata = 4) # algorithm will produce 4 strata

#--- apply stratified sampling algorithm ---#
existing <- sample_strat(sraster = sraster, # use mraster as input for sampling
                         nSamp = 200, # request 200 samples be taken
                         mindist = 100) # define that samples must be 100 m apart


## ----warning=F,message=F------------------------------------------------------
#--- quantile sraster ---#
quantiles <- strat_quantiles(mraster = mraster$zq90,
                             nStrata = 8)

#--- random samples ---#
srs <- sample_srs(raster = sraster,
                  nSamp = 50)

#--- calculate representation ---#
calculate_representation(sraster = quantiles, 
                         existing = srs, 
                         plot = TRUE)

## ----warning=F,message=F------------------------------------------------------
calculate_representation(sraster = sraster,
                         existing = existing,
                         plot = TRUE)

## ----warning=F,message=F------------------------------------------------------
calculate_distance(raster = sraster, # input
                   access = access, # define access road network
                   plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
calculate_pcomp(mraster = mraster, # input
                nComp = 3, # number of components to output
                plot = TRUE, # plot
                details = TRUE) # details about the principal component analysis appended

## ---- warning = FALSE---------------------------------------------------------
#--- determine sample size based on relative standard error (rse) of 1% ---#
calculate_sampsize(mraster = mraster,
                   rse = 0.01)


## ---- warning = FALSE---------------------------------------------------------
#--- change default threshold sequence values ---# 
#--- if increment and rse are not divisible the closest value will be taken ---#
p <- calculate_sampsize(mraster = mraster,
                   rse = 0.025,
                   start = 0.01,
                   end = 0.08,
                   increment = 0.01,
                   plot = TRUE)

p

## ----warning=F,message=F------------------------------------------------------
#--- perform grid sampling ---#
calculate_allocation(sraster = sraster, 
                     nSamp = 200)

## ----warning=F,message=F------------------------------------------------------
#--- calculate existing samples to include ---#
e.sr <- extract_strata(sraster = sraster, 
                       existing = existing)

calculate_allocation(sraster = sraster, 
                     nSamp = 200, 
                     existing = e.sr)

## ---- warning=F,message=F-----------------------------------------------------
calculate_allocation(sraster = sraster, # stratified raster
                     nSamp = 200, # desired sample number
                     existing = e.sr, #existing samples
                     allocation = "optim", # optimal allocation
                     mraster = mraster$zq90, # metric raster
                     force = TRUE) # force nSamp number


## -----------------------------------------------------------------------------
calculate_allocation(sraster = sraster, # stratified raster
                     nSamp = 20, # desired sample number
                     allocation = "equal") # optimal allocation

## -----------------------------------------------------------------------------
weights <- c(0.2, 0.2, 0.2, 0.4)

calculate_allocation(sraster = sraster, # stratified raster
                     nSamp = 20, # desired sample number
                     allocation = "manual", # manual allocation
                     weights = weights) # weights adding to 1


## ----warning=F,message=F, eval = FALSE----------------------------------------
#  calculate_coobs(mraster = mraster, # input
#                  existing = existing, # existing samples
#                  cores = 4, # parallel cores to use
#                  details = TRUE, # provide details from algorithm output
#                  plot = TRUE) # plot

## ----warning=F,message=F, eval = FALSE----------------------------------------
#  #--- by default all statistical data are calculated ---#
#  calculate_pop(mraster = mraster) # input

## ----warning=F,message=F, eval = FALSE----------------------------------------
#  #--- statistical analyses can be chosen by setting their parameter to `FALSE` ---#
#  mat <- calculate_pop(mraster = mraster, # input
#                       nQuant = 10) # desired number of quantiles
#  
#  #--- use matrix output within sample ahels algorithm ---#
#  sample_ahels(mraster = mraster,
#               existing = existing,
#               nQuant = 10,
#               nSamp = 50,
#               matCov = mat)

## ----warning=F,message=F, eval = FALSE----------------------------------------
#  #--- calculate lhsPop details ---#
#  poplhs <- calculate_pop(mraster = mr)
#  
#  calculate_lhsOpt(popLHS = poplhs)

## ----warning=F,message=F, eval = FALSE----------------------------------------
#  calculate_lhsOpt(popLHS = poplhs,
#                   PCA = FALSE,
#                   iter = 200)

