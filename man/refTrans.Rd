\name{refTrans}

\alias{refTrans}

\title{
Transform reference scale of data
}

\description{
This function conducts transformations to convert isotope measurements between reference scales.
}

\usage{
refTrans(samples, marker = "d2H", ref_scale = "VSMOW_H", niter = 5000)
}

\arguments{
  \item{samples}{
  data.frame. Must include a field with data to be transformed, analytical reproducibility of sample data (1 standard deviation), and original reference scale for calibration of data. These fields must be named \emph{marker}, \emph{marker.sd}, and \emph{marker_cal}, respectively, where \emph{marker} is \dQuote{d2H} or \dQuote{d18O}. Values for the cal field should correspond to \emph{Calibration} codes found in \code{\link[assignR]{stds}} tables \emph{hstds} and \emph{ostds}.
}
  \item{marker}{
character string. Column name for isotopic data to be extracted, either \dQuote{d2H} or \dQuote{d18O}.
}
  \item{ref_scale}{
character string. Text identifier for reference scale to which all isotope values will be transformed. See \code{\link[assignR]{stds}}.}
  \item{niter}{
integer. Number of random samples used to propagate uncertainty in calibration hierarchy transformations.
}
}

\value{
Returns an object of class \dQuote{refTrans}.
  \item{data}{
data.frame. Formatted identically to input object \emph{samples}, with values for the data and data uncertainty fields replaced with transformed values.}
  \item{chains}{
list. Each item is a character string vector containing the hierarchy of calibrations used in the transformation for a set of selected samples. See \code{\link[assignR]{stds}}.
}
}

\references{
Magozzi, S. et al. (in press) Calibration chain transformation to improve the comparability of organic hydrogen and oxygen isotope data. \emph{Methods in Ecology and Evolution}
}

\examples{
# Some fake sample data
s = data.frame("d2H" = seq(-100, -10, by=10), "d2H.sd" = rep(2), "d2H_cal" = rep("OldUT_H_1"))

# Transform to VSMOW-SLAP scale using default arguments
d1 = refTrans(s)

# Transformed values
d1$data$d2H

# error - target scale not valid for marker
\dontrun{d2 = refTrans(s, ref_scale = "VSMOW_O")
}
}
