id = letters[1:5]
set.seed(123)
d2H = rnorm(5, -110, 8)
d2H.sd = runif(5, 1.5, 2.5)
d2H_cal = rep("UT_H_1", 5)
un1 = data.frame(id, d2H, d2H.sd, d2H_cal)
un2 = data.frame(id, "sam" = d2H, d2H.sd, d2H_cal)
un3 = data.frame(id, d2H, "sam" = d2H.sd, d2H_cal)
un4 = data.frame(id, d2H, d2H.sd, "sam" = d2H_cal)
un5 = data.frame(id, d2H, d2H.sd, "d2H_cal" = rep("sam", 5))
un6 = data.frame(id, "d18O" = d2H, "d18O.sd" = d2H.sd, 
                 "d18O_cal" = rep("UT_O_4", 5))
r = refTrans(un1, niter = 100)

test_that("refTrans can correctly transform data",{
  expect_s3_class(r, "refTrans")
  expect_equal(length(r[[2]]), 1)
  expect_error(refTrans(un2))
  expect_error(refTrans(un3))
  expect_error(refTrans(un4))
  expect_error(suppressWarnings(refTrans(un5)))
  expect_error(refTrans(un1, marker = "d18O"))
  expect_s3_class(refTrans(un6, marker = "d18O", ref_scale = "VSMOW_O",
                           niter = 100), "refTrans")
})
  