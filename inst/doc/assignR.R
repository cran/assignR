## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)

## ----install, message=FALSE, warning=FALSE, results="hide"--------------------
library(assignR)
library(raster)
library(sp)

## ----boundary-----------------------------------------------------------------
data("naMap")
plot(naMap)

## ----isoscape, fig.width=7, fig.asp=0.45--------------------------------------
data("d2h_lrNA")
plot(d2h_lrNA)

## ----knownOrig_names----------------------------------------------------------
data("knownOrig")
names(knownOrig$sites)
names(knownOrig$samples)
names(knownOrig$sources)

## ----knownOrig_sites, fig.width=5, fig.asp=0.8--------------------------------
plot(assignR:::wrld_simpl)
points(knownOrig$sites, col = "red")

## ----knownOrig_taxa-----------------------------------------------------------
unique(knownOrig$samples$Taxon)

## ----birdData, fig.width=5, fig.asp=0.8---------------------------------------
d = subOrigData(taxon = "Lanius ludovicianus", mask = naMap)

## ----birdChains---------------------------------------------------------------
d$chains

## ----birdSources--------------------------------------------------------------
d$sources[,1:3]

## ----birdNoTrans, fig.width=5, fig.asp=0.8------------------------------------
d = subOrigData(taxon = "Lanius ludovicianus", mask = naMap, ref_scale = NULL)
d$sources$H_std_scale

## ----calRaster, fig.width=5, fig.asp=0.8, out.width='45%'---------------------
r = calRaster(known = d, isoscape = d2h_lrNA, mask = naMap)

## ----samples------------------------------------------------------------------
id = letters[1:5]
set.seed(123)
d2H = rnorm(5, -110, 8)
d2H.sd = runif(5, 1.5, 2.5)
d2H_cal = rep("UT_H_1", 5)
un = data.frame(id, d2H, d2H.sd, d2H_cal)
print(un)

## ----refTrans-----------------------------------------------------------------
un = refTrans(un, ref_scale = "OldEC.1_H_1")
print(un)

## ----pdRaster, fig.width=5, fig.asp=0.8, out.width='45%'----------------------
asn = pdRaster(r, unknown = un$data)

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

## ----qtlRaster1, fig.width=5, fig.asp=0.8, out.width='45%'--------------------
qtlRaster(asn, threshold = 0.1)

## ----qtlRaster2, fig.width=5, fig.asp=0.8, out.width='45%'--------------------
qtlRaster(asn, threshold = 0.8, thresholdType = "prob")

## ----jointP, fig.width=5, fig.asp=0.8-----------------------------------------
jointP(asn)

## ----unionP, fig.width=5, fig.asp=0.8-----------------------------------------
up = unionP(asn)

## ----qtlRaster3, fig.width=5, fig.asp=0.8-------------------------------------
qtlRaster(up, threshold = 0.1)

## ----QA1, warning=FALSE-------------------------------------------------------
qa1 = QA(d, d2h_lrNA, valiStation = 8, valiTime = 4, by = 5, mask = naMap, name = "normal")

## ----plot.QA1, fig.width=4, fig.asp=1, out.width='45%'------------------------
plot(qa1)

## ----modraster, fig.width=5, fig.asp=0.8--------------------------------------
dv = getValues(d2h_lrNA[[1]])
dv = dv + rnorm(length(dv), 0, 15)
d2h_fuzzy = setValues(d2h_lrNA[[1]], dv)
plot(d2h_fuzzy)

## ----QA2, warning=FALSE-------------------------------------------------------
d2h_fuzzy = brick(d2h_fuzzy, d2h_lrNA[[2]])
qa2 = QA(d, d2h_fuzzy, valiStation = 8, valiTime = 4, by = 5, mask = naMap, name = "fuzzy")

## ----plot.QA2, fig.width=4, fig.asp=1, out.width='45%'------------------------
plot(qa1, qa2)

