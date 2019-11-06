## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, message=FALSE, warning=FALSE, results="hide"---------------
library(assignR)

## ----boundary------------------------------------------------------------
data("naMap")
plot(naMap)

## ----isoscape------------------------------------------------------------
data("d2h_world")
d2h_coarse = aggregate(d2h_world, 9)
plot(d2h_coarse)

## ----knownOrig_names-----------------------------------------------------
data("knownOrig")
names(knownOrig)

## ----knownOrig_taxa------------------------------------------------------
unique(knownOrig$Taxon)

## ----birdData------------------------------------------------------------
d = subOrigData(taxon = "Lanius ludovicianus", reference = "Hobson et al. 2012", mask = naMap)

## ----calRaster-----------------------------------------------------------
r = calRaster(known = d, isoscape = d2h_coarse, mask = naMap)

## ----samples-------------------------------------------------------------
id = letters[1:5]
set.seed(123)
d2H = rnorm(5, -110, 8)
un = data.frame(id, d2H)
print(un)

## ----pdRaster------------------------------------------------------------
asn = pdRaster(r, unknown = un)

## ----sums----------------------------------------------------------------
cellStats(asn[[1]], 'sum')

## ----polygons------------------------------------------------------------
data("plover_range_BreedingSeason")
data("plover_range_NonBreedingSeason")
plot(naMap)
lines(plover_range_BreedingSeason, col = c("red"))
lines(plover_range_NonBreedingSeason, col = c("blue"))

## ----oddsRatio1----------------------------------------------------------
p12 = rbind(plover_range_BreedingSeason, plover_range_NonBreedingSeason)
oddsRatio(asn, p12)

## ----oddsRatio2----------------------------------------------------------
pp1 = c(-108,42)
pp2 = c(-103,25)
pp12 = SpatialPoints(coords = rbind(pp1,pp2), proj4string=crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
oddsRatio(asn, pp12)

## ----qtlRaster1----------------------------------------------------------
qtlRaster(asn, threshold = 0.1)

## ----qtlRaster2----------------------------------------------------------
qtlRaster(asn, threshold = 0.8, thresholdType = "prob")

## ----jointP--------------------------------------------------------------
jointP(asn)

## ----unionP--------------------------------------------------------------
up = unionP(asn)

## ----qtlRaster3----------------------------------------------------------
qtlRaster(up, threshold = 0.1)

## ----QA1, warning=FALSE--------------------------------------------------
qa1 = QA(d2h_coarse, d, valiStation = 8, valiTime = 4, mask = naMap, name = "normal")

## ----plot.QA1------------------------------------------------------------
plot(qa1)

## ----modraster-----------------------------------------------------------
dv = getValues(d2h_coarse[[1]])
dv = dv + rnorm(length(dv), 0, 15)
d2h_fuzzy = setValues(d2h_coarse[[1]], dv)
plot(d2h_fuzzy)

## ----QA2, warning=FALSE--------------------------------------------------
d2h_fuzzy = brick(d2h_fuzzy, d2h_coarse[[2]])
qa2 = QA(d2h_fuzzy, d, valiStation = 8, valiTime = 4, mask = naMap, name = "fuzzy")

## ----plot.QA2------------------------------------------------------------
plot(qa1, qa2)

