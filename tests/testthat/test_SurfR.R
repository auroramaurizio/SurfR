test_that("Gene2SProtein -package core function- tests", {

  # 1: Check the output for a specific SP gene
  result1 <- Gene2SProtein(c("EPCAM"), output_tsv = FALSE)
  expect_equal(result1$Surfaceome.Label, "surface")

  # 2: Check the output for an unknown gene, not SP coding
  expect_warning(Gene2SProtein(c("InventedGene"), output_tsv = FALSE))

})

test_that("metaRNAseq", {
  data(ind_deg)
  # 1: if test_statistic is invnorm check that ind_deg and nrep and have the
  # same number of elements
  expect_no_error(metaRNAseq(ind_deg, test_statistic = "invnorm", BHth = 0.05, nrep = c(2,2)))

  # 2: if test_statistic is invnorm nrep must be specified
  expect_error(metaRNAseq(ind_deg, test_statistic = "invnorm", BHth = 0.05))

  # 3: Check that input ind_deg is a list of at least two data.frames.
  # stop when not.
  expect_error(metaRNAseq(ind_deg$DEG2_df, test_statistic = "fishercomb", BHth = 0.05))

  # 4: Check that input ind_deg is a list of at least two data.frames.
  # continue when it is.
  expect_no_error(metaRNAseq(ind_deg, test_statistic = "fishercomb", BHth = 0.05))
})
