## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)

## ----load, message=FALSE, warning=FALSE, results="hide"-----------------------
library(assignR)
library(terra)

## ----boundary-----------------------------------------------------------------
plot(naMap)

## ----isoscape, fig.width=7, fig.asp=0.45--------------------------------------
plot(d2h_lrNA)

## ----knownOrig_names----------------------------------------------------------
names(knownOrig$sites)
names(knownOrig$samples)
names(knownOrig$sources)

## ----knownOrig_sites, fig.width=6, fig.asp=0.6--------------------------------
plot(wrld_simpl)
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

## ----calRaster, fig.width=6, fig.asp=0.8, out.width='90%'---------------------
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

## ----pdRaster, fig.width=6, fig.asp=0.6, out.width='95%'----------------------
Ll_prob = pdRaster(d2h_Ll, Ll_un)

## ----sums---------------------------------------------------------------------
global(Ll_prob[[1]], 'sum', na.rm = TRUE)

## ----Dp, fig.width=5, fig.asp=0.8, out.width='45%'----------------------------
Dp_d = subOrigData(taxon = "Danaus plexippus")
d2h_Dp = calRaster(Dp_d, d2h_lrNA)

## ----srIso, fig.width=5, fig.asp=0.8, out.width='45%'-------------------------
plot(sr_MI$weathered.mean)
crs(sr_MI, describe = TRUE)
crs(d2h_Dp$isoscape.rescale, describe = TRUE)

## ----isoStack-----------------------------------------------------------------
Dp_multi = isoStack(d2h_Dp, sr_MI)
lapply(Dp_multi, crs, describe = TRUE)

## ----Dp_unknown---------------------------------------------------------------
Dp_unk = data.frame("ID" = c("A", "B"), "d2H" = c(-86, -96), "Sr" = c(0.7089, 0.7375))

## ----Dp_Honly, fig.width=5, fig.asp=0.6, out.width='85%'----------------------
Dp_pd_Honly = pdRaster(Dp_multi[[1]], Dp_unk[,-3])

## ----Dp_multi, fig.width=5, fig.asp=0.6, out.width='85%'----------------------
Dp_pd_multi = pdRaster(Dp_multi, Dp_unk)

## ----polygons-----------------------------------------------------------------
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
pp12 = vect(rbind(pp1,pp2))
crs(pp12) = crs(naMap)
oddsRatio(Ll_prob, pp12)

## ----wDist1, fig.width=5, fig.asp=0.8, out.width='45%'------------------------
# View the data
plot(Ll_prob[[1]], main = names(Ll_prob)[1])
points(pp12[1], cex = 2)
plot(Ll_prob[[2]], main = names(Ll_prob)[2])
points(pp12[2], cex = 2)

## ----wDist2, fig.width=5, fig.asp=0.8, out.width='45%'------------------------
wd = wDist(Ll_prob[[1:2]], pp12)
c(wd)[c(1,2,4,6,8,10,12,14,16)] #only showing select columns for formatting!
plot(wd)

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
dv = values(d2h_lrNA[[1]])
dv = dv + rnorm(length(dv), 0, 15)
d2h_fuzzy = setValues(d2h_lrNA[[1]], dv)
plot(d2h_fuzzy)

## ----QA2, warning=FALSE, results='hide'---------------------------------------
d2h_fuzzy = c(d2h_fuzzy, d2h_lrNA[[2]])
qa2 = QA(Ll_d, d2h_fuzzy, valiStation = 8, valiTime = 4, by = 5, mask = naMap, name = "fuzzy")

## ----plot.QA2, fig.width=4, fig.asp=1, out.width='45%'------------------------
plot(qa1, qa2)

