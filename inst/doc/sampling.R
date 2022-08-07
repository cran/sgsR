## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=F,message=F,echo=FALSE-------------------------------------------
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
#--- perform simple random sampling ---#
sample_srs(raster = sraster, # input sraster
           nSamp = 200, # number of desired samples
           plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
sample_srs(raster = mraster, # input mraster
           nSamp = 200, # number of desired samples
           access = access, # define access road network
           mindist = 200, # minimum distance samples must be apart from one another
           buff_inner = 50, # inner buffer - no samples within this distance from road
           buff_outer = 200, # outer buffer - no samples further than this distance from road
           plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
#--- perform grid sampling ---#
sample_systematic(raster = sraster, # input sraster
                  cellsize = 1000, # grid distance
                  plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
#--- perform grid sampling ---#
sample_systematic(raster = sraster, # input sraster
                  cellsize = 500, # grid distance
                  square = FALSE, # hexagonal tessellation
                  location = "random", # random sample within tessellation
                  plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
sample_systematic(raster = sraster, # input sraster
            cellsize = 500, # grid distance
            access = access, # define access road network
            buff_outer = 200, # outer buffer - no samples further than this distance from road
            square = FALSE, # hexagonal tessellation
            location = "corners", # take corners instead of centers
            plot = TRUE)

## ----warning=F,message=F------------------------------------------------------
#--- perform stratified sampling random sampling ---#
sample_strat(sraster = sraster, # input sraster
             nSamp = 200, # desired sample number
             plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
#--- extract strata values to existing samples ---#              
e.sr <- extract_strata(sraster = sraster, # input sraster
                       existing = existing) # existing samples to add strata value to


## ----warning=F,message=F------------------------------------------------------
sample_strat(sraster = sraster, # input sraster
             nSamp = 200, # desired sample number
             access = access, # define access road network
             existing = e.sr, # existing samples with strata values
             mindist = 200, # minimum distance samples must be apart from one another
             buff_inner = 50, # inner buffer - no samples within this distance from road
             buff_outer = 200, # outer buffer - no samples further than this distance from road
             plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
sample_strat(sraster = sraster, # input
             nSamp = 200, # desired sample number
             access = access, # define access road network
             existing = e.sr, # existing samples with strata values
             include = TRUE, # include existing plots in nSamp total
             buff_outer = 200, # outer buffer - no samples further than this distance from road
             plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
#--- perform stratified sampling random sampling ---#
sample_strat(sraster = sraster, # input sraster
             method = "random", #stratified random sampling
             nSamp = 200, # desired sample number
             plot = TRUE) # plot

## -----------------------------------------------------------------------------
#--- perform simple random sampling ---#
sample_nc(mraster = mraster, # input
          nSamp = 25, # desired sample number
          plot = TRUE)

## -----------------------------------------------------------------------------
#--- perform simple random sampling ---#
samples <- sample_nc(mraster = mraster, # input
                    k = 2, # number of nearest neighbours to take for each kmeans center
                    nSamp = 25, # desired sample number
                    plot = TRUE)

#--- total samples = nSamp * k (25 * 2) = 50 ---#
nrow(samples)

## -----------------------------------------------------------------------------
#--- perform simple random sampling with details ---#
details <- sample_nc(mraster = mraster, # input
                     nSamp = 25, # desired sample number
                     details = TRUE)

#--- plot ggplot output ---#

details$kplot

## ----eval = FALSE-------------------------------------------------------------
#  sample_clhs(mraster = mraster, # input
#              nSamp = 200, # desired sample number
#              plot = TRUE, # plot
#              iter = 100) # number of iterations

## ----warning=F,message=F,echo=F,results = FALSE-------------------------------
sample_clhs(mraster = mraster, # input
            nSamp = 200, # desired sample number
            plot = TRUE, # plot 
            iter = 100) # number of iterations

## ----eval = FALSE-------------------------------------------------------------
#  sample_clhs(mraster = mraster, # input
#              nSamp = 300, # desired sample number
#              iter = 100, # number of iterations
#              existing = existing, # existing samples
#              access = access, # define access road network
#              buff_inner = 100, # inner buffer - no samples within this distance from road
#              buff_outer = 300, # outer buffer - no samples further than this distance from road
#              plot = TRUE) # plot

## ----warning=F,message=F,echo=F,results = FALSE-------------------------------
sample_clhs(mraster = mraster, # input
            nSamp = 300, # desired sample number
            iter = 100, # number of iterations
            existing = existing, # existing samples
            access = access, # define access road network
            buff_inner = 100, # inner buffer - no samples within this distance from road
            buff_outer = 300, # outer buffer - no samples further than this distance from road
            plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
#--- cost constrained examples ---#
#--- calculate distance to access layer for each pixel in mr ---#
mr.c <- calculate_distance(raster = mraster, # input
                           access = access,
                           plot = TRUE) # define access road network


## ----eval=F-------------------------------------------------------------------
#  sample_clhs(mraster = mr.c, # input
#              nSamp = 250, # desired sample number
#              iter = 100, # number of iterations
#              cost = "dist2access", # cost parameter - name defined in calculate_distance()
#              plot = TRUE) # plot

## ----warning=F,message=F,echo=F,results = FALSE-------------------------------
sample_clhs(mraster = mr.c, # input
            nSamp = 250, # desired sample number
            iter = 100, # number of iterations
            cost = "dist2access", # cost parameter - name defined in calculate_distance()
            plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
sample_balanced(mraster = mraster, # input
                nSamp = 200, # desired sample number
                plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
sample_balanced(mraster = mraster, # input
                nSamp = 100, # desired sample number
                algorithm = "lcube", # algorithm type
                access = access, # define access road network
                buff_inner = 50, # inner buffer - no samples within this distance from road
                buff_outer = 200) # outer buffer - no samples further than this distance from road

## ----eval = FALSE-------------------------------------------------------------
#  sample_ahels(mraster = mraster,
#               existing = existing, # existing samples
#               plot = TRUE) # plot

## ----warning=F,message=F,echo=FALSE, results = FALSE--------------------------
s <- sample_ahels(mraster = mraster, 
             existing = existing, # existing samples
             plot = TRUE) # plot

## ----echo=FALSE---------------------------------------------------------------
s

## ----eval = FALSE-------------------------------------------------------------
#  sample_ahels(mraster = mraster,
#               existing = existing, # existing samples
#               nQuant = 20, # define 20 quantiles
#               nSamp = 300) # total samples desired

## ----warning=F,message=F,echo=FALSE, results = FALSE--------------------------
s <- sample_ahels(mraster = mraster, 
             existing = existing, # existing samples
             nQuant = 20, # define 20 quantiles
             nSamp = 300, # total samples desired
             plot = TRUE) # plot


## ----echo=FALSE---------------------------------------------------------------
s

