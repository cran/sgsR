% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_strata.R
\name{extract_strata}
\alias{extract_strata}
\title{Extract strata}
\usage{
extract_strata(
  sraster,
  existing,
  quiet = FALSE,
  data.frame = FALSE,
  filename = NULL,
  overwrite = FALSE
)
}
\arguments{
\item{sraster}{spatRaster. Stratification raster.}

\item{existing}{sf 'POINT'.  Existing plot network.}

\item{quiet}{Logical. If \code{TRUE} the user will not get messages
about samples with \code{NA} values.}

\item{data.frame}{Logical. Output as data.frame if \code{TRUE}}

\item{filename}{Character. Path to write output samples.}

\item{overwrite}{Logical. Choice to overwrite existing \code{filename} if it exists.}
}
\value{
An sf or data.frame object of samples with strata attribute.
}
\description{
Extract stratum values to existing samples
}
\note{
If \code{data.frame = TRUE} output will be written using \code{\link[utils]{write.table}}
}
\examples{
#--- Load sraster ---#
r <- system.file("extdata", "sraster.tif", package = "sgsR")
sr <- terra::rast(r)

#--- load existing samples ---#
e <- system.file("extdata", "existing.shp", package = "sgsR")
e <- sf::st_read(e)

extract_strata(
  sraster = sr,
  existing = e
)

}
\seealso{
Other extract functions: 
\code{\link{extract_metrics}()}
}
\author{
Tristan R.H. Goodbody
}
\concept{extract functions}
