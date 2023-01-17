## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=F,message=F,echo=FALSE-------------------------------------------
library(sgsR)
library(terra)
library(dplyr)

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

#--- algorithm table ---#

a <- c("`sample_srs()`","`sample_systematic()`","`sample_strat()`","`sample_nc()`","`sample_clhs()`","`sample_balanced()`","`sample_ahels()`","`sample_existing()`")

d <- c("Simple random", "Systematic", "Stratified", "Nearest centroid", "Conditioned Latin hypercube", "Balanced sampling", "Adapted hypercube evaluation of a legacy sample", "Sub-sampling an existing sample")

s <- c("", "", "[Queinnec, White, & Coops (2021)](https://www.sciencedirect.com/science/article/pii/S0034425721002303)","[Melville & Stone (2016)](https://doi.org/10.1080/00049158.2016.1218265)", "[Minasny & McBratney (2006)](https://www.sciencedirect.com/science/article/pii/S009830040500292X?via%3Dihub)","[Grafström, A. Lisic, J (2018)](http://www.antongrafstrom.se/balancedsampling/)","[Malone, Minasny, & Brungard (2019)](https://peerj.com/articles/6451/)","")

urls <- c("#srs","#systematic","#sstrat","#nc","#clhs","#balanced","#ahels","#samp-existing")

df <- data.frame(Algorithm = a, Description = d, Reference = s)

df$Algorithm <- paste0("[", df$Algorithm, "](", urls, ")")


## ----echo=FALSE---------------------------------------------------------------
knitr::kable(df, align = 'c')

## ----warning=F,message=F------------------------------------------------------
#--- perform simple random sampling ---#
sample_srs(raster = sraster, # input sraster
           nSamp = 200, # number of desired sample units
           plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
sample_srs(raster = mraster, # input mraster
           nSamp = 200, # number of desired sample units
           access = access, # define access road network
           mindist = 200, # minimum distance sample units must be apart from one another
           buff_inner = 50, # inner buffer - no sample units within this distance from road
           buff_outer = 200, # outer buffer - no sample units further than this distance from road
           plot = TRUE) # plot

## ----warning=F,message=F,eval = FALSE-----------------------------------------
#  #--- perform grid sampling ---#
#  sample_systematic(raster = sraster, # input sraster
#                    cellsize = 1000, # grid distance
#                    plot = TRUE) # plot

## ----warning=F,message=F,eval = FALSE-----------------------------------------
#  #--- perform grid sampling ---#
#  sample_systematic(raster = sraster, # input sraster
#                    cellsize = 500, # grid distance
#                    square = FALSE, # hexagonal tessellation
#                    location = "random", # randomly sample within tessellation
#                    plot = TRUE) # plot

## ----warning=F,message=F,eval = FALSE-----------------------------------------
#  sample_systematic(raster = sraster, # input sraster
#              cellsize = 500, # grid distance
#              access = access, # define access road network
#              buff_outer = 200, # outer buffer - no sample units further than this distance from road
#              square = FALSE, # hexagonal tessellation
#              location = "corners", # take corners instead of centers
#              plot = TRUE)

## ----warning=F,message=F------------------------------------------------------
#--- perform stratified sampling random sampling ---#
sample_strat(sraster = sraster, # input sraster
             nSamp = 200, # desired sample size
             plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
#--- extract strata values to existing samples ---#              
e.sr <- extract_strata(sraster = sraster, # input sraster
                       existing = existing) # existing samples to add strata value to


## ----warning=F,message=F------------------------------------------------------
sample_strat(sraster = sraster, # input sraster
             nSamp = 200, # desired sample size
             access = access, # define access road network
             existing = e.sr, # existing sample with strata values
             mindist = 200, # minimum distance sample units must be apart from one another
             buff_inner = 50, # inner buffer - no sample units within this distance from road
             buff_outer = 200, # outer buffer - no sample units further than this distance from road
             plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
sample_strat(sraster = sraster, # input
             nSamp = 200, # desired sample size
             access = access, # define access road network
             existing = e.sr, # existing samples with strata values
             include = TRUE, # include existing sample in nSamp total
             buff_outer = 200, # outer buffer - no samples further than this distance from road
             plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
#--- perform stratified sampling random sampling ---#
sample_strat(sraster = sraster, # input sraster
             method = "random", #stratified random sampling
             nSamp = 200, # desired sample size
             plot = TRUE) # plot

## -----------------------------------------------------------------------------
#--- perform simple random sampling ---#
sample_nc(mraster = mraster, # input
          nSamp = 25, # desired sample size
          plot = TRUE)

## -----------------------------------------------------------------------------
#--- perform simple random sampling ---#
samples <- sample_nc(mraster = mraster, # input
                    k = 2, # number of nearest neighbours to take for each kmeans center
                    nSamp = 25, # desired sample size
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
#              nSamp = 200, # desired sample size
#              plot = TRUE, # plot
#              iter = 100) # number of iterations

## ----warning=F,message=F,echo=F,results = FALSE-------------------------------
sample_clhs(mraster = mraster, # input
            nSamp = 200, # desired sample size
            plot = TRUE, # plot 
            iter = 100) # number of iterations

## ----warning=F,message=F------------------------------------------------------
#--- cost constrained examples ---#
#--- calculate distance to access layer for each pixel in mr ---#
mr.c <- calculate_distance(raster = mraster, # input
                           access = access, # define access road network
                           plot = TRUE) # plot


## ----eval=F-------------------------------------------------------------------
#  sample_clhs(mraster = mr.c, # input
#              nSamp = 250, # desired sample size
#              iter = 100, # number of iterations
#              cost = "dist2access", # cost parameter - name defined in calculate_distance()
#              plot = TRUE) # plot

## ----warning=F,message=F,echo=F,results = FALSE-------------------------------
sample_clhs(mraster = mr.c, # input
            nSamp = 250, # desired sample size
            iter = 100, # number of iterations
            cost = "dist2access", # cost parameter - name defined in calculate_distance()
            plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
sample_balanced(mraster = mraster, # input
                nSamp = 200, # desired sample size
                plot = TRUE) # plot

## ----warning=F,message=F------------------------------------------------------
sample_balanced(mraster = mraster, # input
                nSamp = 100, # desired sample size
                algorithm = "lcube", # algorithm type
                access = access, # define access road network
                buff_inner = 50, # inner buffer - no sample units within this distance from road
                buff_outer = 200) # outer buffer - no sample units further than this distance from road

## ----eval = FALSE-------------------------------------------------------------
#  #--- remove `type` variable from existing  - causes plotting issues ---#
#  
#  existing <- existing %>% select(-type)
#  
#  sample_ahels(mraster = mraster,
#               existing = existing, # existing sample
#               plot = TRUE) # plot

## ----warning=F,message=F,echo=FALSE, results = FALSE--------------------------
s <- sample_ahels(mraster = mraster, 
             existing = existing) # existing samples

## ----echo=FALSE---------------------------------------------------------------
s

## ----eval = FALSE-------------------------------------------------------------
#  sample_ahels(mraster = mraster,
#               existing = existing, # existing sample
#               nQuant = 20, # define 20 quantiles
#               nSamp = 300) # desired sample size

## ----warning=F,message=F,echo=FALSE, results = FALSE--------------------------
s <- sample_ahels(mraster = mraster, 
             existing = existing, # existing sample
             nQuant = 20, # define 20 quantiles
             nSamp = 300) # plot


## ----echo=FALSE---------------------------------------------------------------
s

## ----warning=F,message=F------------------------------------------------------
#--- generate existing samples and extract metrics ---#
existing <- sample_srs(raster = mraster, nSamp = 900, plot = TRUE)


## ----warning=F,message=F------------------------------------------------------
#--- sub sample using ---#
e <- existing %>% 
  extract_metrics(mraster = mraster, existing = .)

sample_existing(existing = e, nSamp = 300, plot = TRUE)


## ----warning=F,message=F------------------------------------------------------
#--- sub sample using ---#
sample_existing(existing = existing, # our existing sample
                nSamp = 300, # desired sample size
                raster = mraster, # include mraster metrics to guide sampling of existing
                plot = TRUE) # plot


## ----warning=F,message=F------------------------------------------------------
#--- create distance from roads metric ---#
dist <- calculate_distance(raster = mraster, access = access)


## ----warning=F,message=F------------------------------------------------------
#--- sub sample using ---#
sample_existing(existing = existing, # our existing sample
                nSamp = 300, # desired sample size
                raster = dist, # include mraster metrics to guide sampling of existing
                cost = 4, # either provide the index (band number) or the name of the cost layer
                plot = TRUE) # plot


## ----warning=F,message=F------------------------------------------------------
#--- ensure access and existing are in the same CRS ---#

sf::st_crs(existing) <- sf::st_crs(access)

#--- sub sample using ---#
sample_existing(existing = existing, # our existing sample
                nSamp = 300, # desired sample size
                raster = dist, # include mraster metrics to guide sampling of existing
                cost = 4, # either provide the index (band number) or the name of the cost layer
                access = access, # roads layer
                buff_inner = 50, # inner buffer - no sample units within this distance from road
                buff_outer = 300, # outer buffer - no sample units further than this distance from road
                plot = TRUE) # plot


