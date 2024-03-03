## ----include = FALSE----------------------------------------------------------
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
sraster <- strat_quantiles(
  mraster = mraster$zq90, # use mraster as input for sampling
  nStrata = 4
) # algorithm will produce 4 strata

#--- apply stratified sampling algorithm ---#
existing <- sample_strat(
  sraster = sraster, # use mraster as input for sampling
  nSamp = 200, # request 200 samples be taken
  mindist = 100
) # define that samples must be 100 m apart

#--- algorithm table ---#

a <- c("`strat_kmeans()`", "`strat_quantiles()`", "`strat_breaks()`", "`strat_poly()`", "`strat_map()`")

d <- c("kmeans", "Quantiles", "User-defined breaks", "Polygons", "Maps (combines) `srasters`")

s <- c("Unsupervised", "Either", "Supervised", "Supervised", "Unsupervised")

urls <- c("#kmeans", "#quantiles", "#breaks", "#poly", "#map")

df <- data.frame(Algorithm = a, Description = d, Approach = s)

df$Algorithm <- paste0("[", df$Algorithm, "](", urls, ")")

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(df, align = "c")

## ----warning=F,message=F------------------------------------------------------
#--- perform stratification using k-means ---#
strat_kmeans(
  mraster = mraster, # input
  nStrata = 5
) # algorithm will produce 4 strata

## ----warning=F,message=F------------------------------------------------------
strat_kmeans(
  mraster = mraster, # input
  nStrata = 10, # algorithm will produce 10 strata
  iter = 1000, # set minimum number of iterations to determine kmeans centers
  algorithm = "MacQueen", # use MacQueen algorithm
  plot = TRUE
) # plot output

## ----warning=F,message=F------------------------------------------------------
#--- perform quantiles stratification ---#
strat_quantiles(
  mraster = mraster$zq90,
  nStrata = 6,
  plot = TRUE
)

#--- vectorized ---#
strat_quantiles(
  mraster = mraster[[1:2]], # two metric layers
  nStrata = list(c(0.2, 0.4, 0.8), 3), # list with two objects - 1 probability breaks, 1 scalar integer
  plot = TRUE, # plot output srasters
  map = TRUE
) # combine stratifications to a mapped output

## ----warning=F,message=F------------------------------------------------------
#--- perform stratification using user-defined breaks ---#

#--- define breaks for metric ---#
br.pz2 <- c(20, 40, 60, 80)

br.pz2

#--- perform stratification using user-defined breaks ---#

#--- define breaks for metric ---#
br.zq90 <- c(3, 5, 11, 18)

br.zq90

## ----warning=F,message=F------------------------------------------------------
#--- stratify on 1 metric only ---#
strat_breaks(
  mraster = mraster$pzabove2, # single raster
  breaks = br.pz2, # single set of breaks
  plot = TRUE
) # plot output

## ----warning=F,message=F------------------------------------------------------
#--- vectorized ---#
strat_breaks(
  mraster = mraster[[1:2]], # two metrics
  breaks = list(br.zq90, br.pz2), # list of two breaks vectors
  map = TRUE, # map final output
  plot = TRUE
) # plot outputs

## -----------------------------------------------------------------------------
#--- load in polygon coverage ---#
poly <- system.file("extdata", "inventory_polygons.shp", package = "sgsR")

fri <- sf::st_read(poly)

#--- specify polygon attribute to stratify ---#

attribute <- "NUTRIENTS"

#--- specify features within attribute & how they should be grouped ---#
#--- as a single vector ---#

features <- c("poor", "rich", "medium")

## -----------------------------------------------------------------------------
#--- stratify polygon coverage ---#

srasterpoly <- strat_poly(
  poly = fri, # input polygon
  attribute = attribute, # attribute to stratify by
  features = features, # features within attribute
  raster = sraster, # raster to define extent and resolution for output
  plot = TRUE
) # plot output

## -----------------------------------------------------------------------------
#--- or as multiple lists ---#
g1 <- "poor"
g2 <- c("rich", "medium")

features <- list(g1, g2)

strat_poly(
  poly = fri,
  attribute = attribute,
  features = features,
  raster = sraster,
  plot = TRUE,
  details = TRUE
)

## -----------------------------------------------------------------------------
#--- stack srasters together ---#

srasters <- c(srasterpoly, sraster)

plot(srasters)

## -----------------------------------------------------------------------------
#--- map srasters ---#
strat_map(
  sraster = srasters, # two layer sraster
  plot = TRUE
)

## -----------------------------------------------------------------------------
strat_map(
  sraster = srasters, # input with 2 sraster layers
  stack = TRUE, # output stacked input (strata_1, strata_2) and output (strata) layers
  details = TRUE, # provide additional details
  plot = TRUE
) # plot output

