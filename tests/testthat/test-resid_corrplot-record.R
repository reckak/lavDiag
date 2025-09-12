test_that("resid_corrplot(record=TRUE) returns recordedplot for single-group", {
  skip_if_no_displaylist()

  HS.model <- '
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
  '
  fit1 <- lavaan::cfa(HS.model, data = lavaan::HolzingerSwineford1939)

  rec <- resid_corrplot(fit1, order = "hclust", hclust.method = "ward.D2", record = TRUE)
  expect_s3_class(rec, "recordedplot")
})

test_that("resid_corrplot(record=TRUE) returns a list of recordedplots for multi-group", {
  skip_if_no_displaylist()

  HS.model <- '
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
  '
  fit2 <- lavaan::cfa(HS.model, data = lavaan::HolzingerSwineford1939, group = "school")

  recs <- resid_corrplot(fit2, order = "hclust", hclust.method = "ward.D", record = TRUE)
  expect_type(recs, "list")

  # number of groups
  ng <- lavaan::lavInspect(fit2, "ngroups")
  expect_length(recs, ng)

  # all elements are recordedplot
  expect_true(all(vapply(recs, inherits, logical(1), what = "recordedplot")))
})
