## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, message=FALSE, warning=FALSE, results="hide"--------------------
library(assignR)
library(raster)
library(sp)

## ----boundary-----------------------------------------------------------------
data("naMap")
plot(naMap)

## ----isoscape, fig.width = 6, fig.asp = 0.5-----------------------------------
data("d2h_lrNA")
plot(d2h_lrNA)

## ----knownOrig_names----------------------------------------------------------
data("knownOrig")
names(knownOrig)

## ----knownOrig_taxa-----------------------------------------------------------
unique(knownOrig$Taxon)

## ----birdData-----------------------------------------------------------------
d = subOrigData(taxon = "Lanius ludovicianus", reference = "Hobson et al. 2012", mask = naMap)

## ----calRaster----------------------------------------------------------------
r = calRaster(known = d, isoscape = d2h_lrNA, mask = naMap)

## ----samples------------------------------------------------------------------
id = letters[1:5]
set.seed(123)
d2H = rnorm(5, -110, 8)
un = data.frame(id, d2H)
print(un)

## ----pdRaster-----------------------------------------------------------------
asn = pdRaster(r, unknown = un)

## ----sums---------------------------------------------------------------------
cellStats(asn[[1]], 'sum')

## ----polygons-----------------------------------------------------------------
data("states")
s1 = states[states$STATE_ABBR == "UT",]
s2 = states[states$STATE_ABBR == "NM",]
plot(naMap)
lines(s1, col = c("red"))
lines(s2, col = c("blue"))

## ----oddsRatio1---------------------------------------------------------------
s12 = rbind(s1, s2)
oddsRatio(asn, s12)

## ----oddsRatio2---------------------------------------------------------------
pp1 = c(-112,40)
pp2 = c(-105,33)
pp12 = SpatialPoints(coords = rbind(pp1,pp2))
proj4string(pp12) = proj4string(naMap)
oddsRatio(asn, pp12)

## ----qtlRaster1---------------------------------------------------------------
qtlRaster(asn, threshold = 0.1)

## ----qtlRaster2---------------------------------------------------------------
qtlRaster(asn, threshold = 0.8, thresholdType = "prob")

## ----jointP-------------------------------------------------------------------
jointP(asn)

## ----unionP-------------------------------------------------------------------
up = unionP(asn)

## ----qtlRaster3---------------------------------------------------------------
qtlRaster(up, threshold = 0.1)

## ----QA1, warning=FALSE-------------------------------------------------------
qa1 = QA(d2h_lrNA, d, valiStation = 8, valiTime = 4, by = 5, mask = naMap, name = "normal")

## ----plot.QA1-----------------------------------------------------------------
plot(qa1)

## ----modraster, fig.width=5, fig.asp = 0.6------------------------------------
dv = getValues(d2h_lrNA[[1]])
dv = dv + rnorm(length(dv), 0, 15)
d2h_fuzzy = setValues(d2h_lrNA[[1]], dv)
plot(d2h_fuzzy)

## ----QA2, warning=FALSE-------------------------------------------------------
d2h_fuzzy = brick(d2h_fuzzy, d2h_lrNA[[2]])
qa2 = QA(d2h_fuzzy, d, valiStation = 8, valiTime = 4, by = 5, mask = naMap, name = "fuzzy")

## ----plot.QA2-----------------------------------------------------------------
plot(qa1, qa2)

