\name{d2h_lrNA}

\alias{d2h_lrNA}

\docType{data}

\title{
data: low resolution North American growing season H isoscape
}

\description{
Interpolated growing season precipitation H isoscape from waterisotopes.org
}

\usage{data("d2h_lrNA")
}

\format{
RasterBrick with two layers. The first layer is the mean prediction and the second is 1 standard deviation
}

\references{
Bowen, G. J. (2018) Gridded maps of the isotopic composition of meteoric waters. http://www.waterisotopes.org.

Bowen, G. J., Wassenaar, L. I. and Hobson, K. A. (2005) Global application of stable hydrogen and oxygen isotopes to wildlife forensics. \emph{Oecologia}, \bold{143}, 337--348.

IAEA/WMO (2018) Global Network of Isotopes in Precipitation. The GNIP Database. https://nucleus.iaea.org/wiser.
}

\examples{
library(raster)

data("d2h_lrNA")
plot(d2h_lrNA)
}

