dmp = capture_output({
  d = subOrigData(taxon = "Homo sapiens", dataset = 10, mask = naMap,
                                   genplot = FALSE)
})

test_that("suOrigData works",{
  expect_equal(class(d), "subOrigData")
  expect_is(d$data, "SpatVector")
  expect_error(subOrigData(taxon = "Turdus philomelos", mask = naMap))
  expect_error(subOrigData(taxon = "Turdus philomelos", marker = "d14C"))
  expect_warning(subOrigData(taxon = "Serin serin", 
                             age_code = c("juvenile", "newborn"),
                             ref_scale = NULL, genplot = FALSE))
  expect_warning(subOrigData(taxon = c("Serin serin", "Vanellus malabaricus"),
                             ref_scale = NULL, genplot = FALSE))
  expect_warning(subOrigData(group = c("Indigenous human", "Badgers"),
                             ref_scale = NULL, genplot = FALSE))
  expect_warning(subOrigData(dataset = c(8, "Ma 2020"),
                             ref_scale = NULL, genplot = FALSE))
  expect_warning(subOrigData(dataset = c(8, 100),
                             ref_scale = NULL, genplot = FALSE))
})

dmp = capture_output({
  d_hasNA = d
  d_hasNA$data$d2H[1] = NA
  d_diffProj = d
  d_diffProj$data = project(d$data, "+init=epsg:28992")
  d_usr_bad = d$data
  d_usr_good = d_usr_bad
  values(d_usr_good) = data.frame(d$data$d2H, d$data$d2H.sd)
  d_noCRS = d
  crs(d_noCRS$data) = ""
  
  d2h_lrNA_noCRS = d2h_lrNA
  crs(d2h_lrNA_noCRS) = ""
  
  mask_diffProj = project(naMap, "+init=epsg:28992")
  
  mask_noCRS = naMap
  crs(mask_noCRS) = ""
  
  tempVals = values(d2h_lrNA)
  tempVals[is.nan(tempVals)] = 9999
  d2h_lrNA_with9999 = setValues(d2h_lrNA, tempVals)
  
  s1 = states[states$STATE_ABBR == "UT",]
  d2h_lrNA_na = mask(d2h_lrNA, s1)
  
  r = calRaster(known = d, isoscape = d2h_lrNA_with9999, NA.value = 9999, 
                interpMethod = 1, genplot = FALSE, mask = naMap)
})

test_that("calRaster works",{
  capture_output({
    expect_is(r, "rescale")
    expect_is(calRaster(known = d_usr_good, isoscape = d2h_lrNA,
                        genplot = FALSE), "rescale")
    expect_output(calRaster(known = d, isoscape = d2h_lrNA, 
                            genplot = FALSE, outDir = tempdir()))
    expect_equal(nlyr(r$isoscape.rescale), 2)
    expect_error(calRaster(known = d$data$d2H, isoscape = d2h_lrNA))
    expect_error(calRaster(known = d, isoscape = d2h_lrNA, 
                           outDir = 2))
    expect_error(calRaster(known = d, isoscape = d2h_lrNA, interpMethod = 3))
    expect_message(calRaster(known = d, isoscape = d2h_lrNA, genplot = 2))
    expect_error(calRaster(known = d, isoscape = d2h_lrNA_noCRS))
    expect_error(calRaster(known = d, isoscape = d2h_lrNA$mean))
    expect_error(calRaster(known = d_usr_bad, isoscape = d2h_lrNA))
    expect_error(calRaster(known = d, isoscape = d2h_lrNA, mask = mask_noCRS))
    expect_error(calRaster(known = d, isoscape = d2h_lrNA, mask = d))
    expect_error(calRaster(known = d_noCRS, isoscape = d2h_lrNA))
    expect_error(calRaster(known = d_hasNA, isoscape = d2h_lrNA, 
                           ignore.NA = FALSE))
    expect_error(calRaster(known = d, isoscape = d2h_lrNA_na, ignore.NA = FALSE))
    expect_message(calRaster(known = d_diffProj, isoscape = d2h_lrNA, 
                             genplot = FALSE))
    expect_message(calRaster(known = d, isoscape = d2h_lrNA, 
                             mask = mask_diffProj, genplot = FALSE))
    expect_warning(calRaster(known = d, isoscape = d2h_lrNA_na, genplot = FALSE))
  })
})

dmp = capture_output({
  id = c("A", "B", "C", "D")
  d2H = c(-110, -90, -105, -102)
  un = data.frame(id,d2H)
  asn = suppressWarnings(pdRaster(r, unknown = un, mask = naMap, genplot = FALSE))
  
  j = jointP(asn)
})

test_that("jointP works",{
  expect_equal(global(j, sum, na.rm = TRUE)[1, 1], 1)
  expect_is(j, "SpatRaster")
  expect_error(jointP(d))
})

dmp = capture_output({
  u = unionP(asn)
})

test_that("unionP works",{
  expect_is(u, "SpatRaster")
  expect_error(unionP(d2H))  
})

dmp = capture_output({
  s1 = states[states$STATE_ABBR == "UT",]
  s2 = states[states$STATE_ABBR == "NM",]
  s12 = rbind(s1, s2)
  o1 = suppressWarnings(oddsRatio(asn, s12))                     
  
  pp1 = c(-112,40)
  pp2 = c(-105,33)
  pp12 = vect(rbind(pp1,pp2), crs = "WGS84")
  o2 = suppressWarnings(oddsRatio(asn, pp12))
  o3 = suppressWarnings(oddsRatio(asn, pp12[1]))
  o4 = suppressWarnings(oddsRatio(asn$A, pp12))
  o5 = suppressWarnings(oddsRatio(asn$A, s12))
  
  s12_diffProj = suppressWarnings(project(s12, "+init=epsg:28992"))
  pp12_diffProj = suppressWarnings(project(pp12, "+init=epsg:28992"))
  
  pp12_noCRS = pp12
  crs(pp12_noCRS) = ""
  s12_noCRS = s12
  crs(s12_noCRS) = ""
  
  pp121 = vect(rbind(pp1, pp2, pp3 = pp1), crs = "WGS84")
})

test_that("oddsRatio works",{
  expect_is(o1, "list")
  expect_is(o2, "list")
  expect_is(o3, "data.frame")
  expect_is(o4, "list")
  expect_is(o5, "list")
  
  expect_error(oddsRatio(naMap,s12))
  expect_error(oddsRatio(asn, data.frame(30.6, 50.5)))
  expect_error(oddsRatio(asn, s12_noCRS))
  expect_error(suppressWarnings(oddsRatio(asn, pp121)))
  expect_error(oddsRatio(asn, s1))
  expect_error(oddsRatio(asn, pp12_noCRS))
  
  expect_message(suppressWarnings(oddsRatio(asn, s12_diffProj)))
  expect_message(suppressWarnings(oddsRatio(asn, pp12_diffProj)))
})
       
dmp = capture_output({
  q1 = qtlRaster(asn, threshold = 0.1, thresholdType = "area", outDir = tempdir())
  q2 = qtlRaster(asn, threshold = 0.1, thresholdType = "prob", genplot = FALSE)
  q3 = qtlRaster(asn, threshold = 0, genplot = FALSE)
})     

test_that("qtlRaster works",{
  expect_is(q1, "SpatRaster")
  expect_equal(nlyr(q1), 4)
  expect_equal(nlyr(q2), 4)
  expect_equal(nlyr(q3), 4)
  expect_error(qtlRaster(asn, threshold = "a"))
  expect_error(qtlRaster(asn, threshold = 10))
  expect_error(qtlRaster(asn, thresholdType = "probability"))
  expect_message(qtlRaster(asn, threshold = 0.1, genplot = "A"))
  expect_error(qtlRaster(asn, threshold = 0.1, outDir = 1))
})

