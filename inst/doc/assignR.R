## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)

## ----load, message=FALSE, warning=FALSE, results="hide"-----------------------
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
Ll_d = subOrigData(taxon = "Lanius ludovicianus", mask = naMap)

## ----birdChains---------------------------------------------------------------
Ll_d$chains

## ----birdSources--------------------------------------------------------------
Ll_d$sources[,1:3]

## ----birdNoTrans, fig.width=5, fig.asp=0.8------------------------------------
Ll_d = subOrigData(taxon = "Lanius ludovicianus", mask = naMap, ref_scale = NULL)
Ll_d$sources$H_cal

## ----calRaster, fig.width=5, fig.asp=0.8, out.width='45%'---------------------
d2h_Ll = calRaster(known = Ll_d, isoscape = d2h_lrNA, mask = naMap)

## ----samples------------------------------------------------------------------
id = letters[1:5]
set.seed(123)
d2H = rnorm(5, -110, 8)
d2H.sd = runif(5, 1.5, 2.5)
d2H_cal = rep("UT_H_1", 5)
Ll_un = data.frame(id, d2H, d2H.sd, d2H_cal)
print(Ll_un)

## ----refTrans-----------------------------------------------------------------
Ll_un = refTrans(Ll_un, ref_scale = "OldEC.1_H_1")
print(Ll_un)

## ----pdRaster, fig.width=5, fig.asp=0.8, out.width='45%'----------------------
Ll_prob = pdRaster(d2h_Ll, unknown = Ll_un)

## ----sums---------------------------------------------------------------------
cellStats(Ll_prob[[1]], 'sum')

## ----Dp, fig.width=5, fig.asp=0.8, out.width='45%'----------------------------
Dp_d = subOrigData(taxon = "Danaus plexippus")
d2h_Dp = calRaster(Dp_d, d2h_lrNA)

## ----srIso, fig.width=5, fig.asp=0.8, out.width='45%'-------------------------
data("sr_MI")
plot(sr_MI$weathered.mean)
proj4string(sr_MI)
proj4string(d2h_Dp$isoscape.rescale)

## ----isoStack-----------------------------------------------------------------
Dp_multi = isoStack(d2h_Dp, sr_MI)
lapply(Dp_multi, proj4string)
lapply(Dp_multi, bbox)

## ----Dp_unknown---------------------------------------------------------------
Dp_unk = data.frame("ID" = c("A", "B"), "d2H" = c(-86, -96), "Sr" = c(0.7089, 0.7375))

## ----Dp_Honly, fig.width=5, fig.asp=0.8, out.width='45%'----------------------
Dp_pd_Honly = pdRaster(Dp_multi[[1]], Dp_unk[,-3])

## ----Dp_multi, fig.width=5, fig.asp=0.8, out.width='45%'----------------------
Dp_pd_multi = pdRaster(Dp_multi, Dp_unk)

## ----polygons-----------------------------------------------------------------
data("states")
s1 = states[states$STATE_ABBR == "UT",]
s2 = states[states$STATE_ABBR == "NM",]
plot(naMap)
plot(s1, col = c("red"), add = TRUE)
plot(s2, col = c("blue"), add = TRUE)

## ----oddsRatio1---------------------------------------------------------------
s12 = rbind(s1, s2)
oddsRatio(Ll_prob, s12)

## ----oddsRatio2---------------------------------------------------------------
pp1 = c(-112,40)
pp2 = c(-105,33)
pp12 = SpatialPoints(coords = rbind(pp1,pp2))
proj4string(pp12) = proj4string(naMap)
oddsRatio(Ll_prob, pp12)

## ----qtlRaster1, fig.width=5, fig.asp=0.8, out.width='45%'--------------------
qtlRaster(Ll_prob, threshold = 0.1)

## ----qtlRaster2, fig.width=5, fig.asp=0.8, out.width='45%'--------------------
qtlRaster(Ll_prob, threshold = 0.8, thresholdType = "prob")

## ----jointP, fig.width=5, fig.asp=0.8-----------------------------------------
jointP(Ll_prob)

## ----unionP, fig.width=5, fig.asp=0.8-----------------------------------------
Ll_up = unionP(Ll_prob)

## ----qtlRaster3, fig.width=5, fig.asp=0.8-------------------------------------
qtlRaster(Ll_up, threshold = 0.1)

## ----QA1, warning=FALSE, results='hide'---------------------------------------
qa1 = QA(Ll_d, d2h_lrNA, valiStation = 8, valiTime = 4, by = 5, mask = naMap, name = "normal")

## ----plot.QA1, fig.width=4, fig.asp=1, out.width='45%'------------------------
plot(qa1)

## ----modraster, fig.width=5, fig.asp=0.8--------------------------------------
dv = getValues(d2h_lrNA[[1]])
dv = dv + rnorm(length(dv), 0, 15)
d2h_fuzzy = setValues(d2h_lrNA[[1]], dv)
plot(d2h_fuzzy)

## ----QA2, warning=FALSE, results='hide'---------------------------------------
d2h_fuzzy = brick(d2h_fuzzy, d2h_lrNA[[2]])
qa2 = QA(Ll_d, d2h_fuzzy, valiStation = 8, valiTime = 4, by = 5, mask = naMap, name = "fuzzy")

## ----plot.QA2, fig.width=4, fig.asp=1, out.width='45%'------------------------
plot(qa1, qa2)

