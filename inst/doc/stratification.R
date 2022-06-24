## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=F,message=F,echo=FALSE,results=FALSE-----------------------------
library(sgsR)
library(terra)
library(sf)

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
#--- perform stratification using k-means ---#
strat_kmeans(mraster = mraster, # input
             nStrata = 5) # algorithm will produce 4 strata


## ----warning=F,message=F------------------------------------------------------
strat_kmeans(mraster = mraster, # input
             nStrata = 10, # algorithm will produce 10 strata
             iter = 1000, # set minimum number of iterations to determine kmeans centers
             algorithm = "MacQueen", # use MacQueen algorithm
             plot = TRUE) # plot output

## ----warning=F,message=F------------------------------------------------------
#--- perform quantiles stratification ---#
strat_quantiles(mraster = mraster$zq90,
                nStrata = 6,
                plot = TRUE)

#--- dual stratification - will produce 12 output strata ---#
strat_quantiles(mraster = mraster$zq90, 
                mraster2 = mraster$zsd,
                nStrata = 3, 
                nStrata2 = 4)

## ----warning=F,message=F------------------------------------------------------
#--- perform stratification using user-defined breaks ---#

#--- define breaks for metric ---#
breaks <- c(seq(0,100,20))

breaks

#--- perform stratification using user-defined breaks ---#

values <- terra::values(mraster$zq90)

#--- define breaks for metric ---#
breaks2 <- c(5,10,15,20,25)

breaks2


## ----warning=F,message=F------------------------------------------------------
#--- stratify on 1 metric only ---#

strat_breaks(mraster = mraster$pzabove2,
             breaks = breaks,
             plot = TRUE)

## ----warning=F,message=F------------------------------------------------------
#--- stratify on 1 metric only ---#

strat_breaks(mraster = mraster$zq90,
             breaks = breaks2,
             plot = TRUE)

## -----------------------------------------------------------------------------
#--- load in polygon coverage ---#
poly <- system.file("extdata", "inventory_polygons.shp", package = "sgsR")

fri <- sf::st_read(poly)

## -----------------------------------------------------------------------------
#--- stratify polygon coverage ---#
#--- specify polygon attribute to stratify ---#

attribute <- "NUTRIENTS"

#--- specify features within attribute & how they should be grouped ---#
#--- as a single vector ---#

features <- c("poor", "rich", "medium")

srasterpoly <- strat_poly(poly = fri, # input polygon
                          attribute = attribute, # attribute to stratify by
                          features = features, # features within attribute
                          raster = sraster, # raster to define extent and resolution for output
                          plot = TRUE) # plot output

## -----------------------------------------------------------------------------
#--- or as multiple lists ---#
g1 <- "poor"
g2 <- c("rich", "medium")

features <- list(g1, g2)

strat_poly(poly = fri,
           attribute = attribute,
           features = features,
           raster = sraster,
           plot = TRUE,
           details = TRUE)

## -----------------------------------------------------------------------------
#--- map srasters ---#
strat_map(sraster = srasterpoly, # strat_poly 3 class stratification
          sraster2 = sraster, # strat_kmeans 4 class stratification
          plot = TRUE)


## -----------------------------------------------------------------------------
strat_map(sraster = srasterpoly, # strat_poly 3 class stratification
          sraster2 = sraster, # strat_poly 3 class stratification
          stack = TRUE, # stack input and output strata into multi layer output raster
          details = TRUE, # provide additional details
          plot = TRUE) # plot output

