data(ind_deg)
data(countData)
data(metadata)
data(enrichedList)

test_that("Gene2SProtein -package core function- tests", {
  # 1: Check the output for a specific SP gene
  result1 <- Gene2SProtein(c("EPCAM"), output_tsv = FALSE)
  expect_equal(result1$Surfaceome.Label, "surface")

  # 2: Check the output for an unknown gene, not SP coding
  expect_warning(Gene2SProtein(c("InventedGene"), output_tsv = FALSE))
})


test_that("Enrichment_barplot", {
  dbs <- c("GO_Cellular_Component_2021")
  expect_no_error(Enrichment_barplot(enrichedList,
                                     enrich.databases = dbs,
                                     p_adj = 0.1, num_term = 2, cond = "UP"))
})

test_that("metaRNAseq", {
  # 1: if test_statistic is invnorm check that ind_deg and nrep and have the
  # same number of elements
  expect_no_error(metaRNAseq(ind_deg,
                             test_statistic = "invnorm",
                             BHth = 0.05,
                             nrep = c(2, 2)))

  # 2: if test_statistic is invnorm nrep must be specified
  expect_error(metaRNAseq(ind_deg,
                          test_statistic = "invnorm",
                          BHth = 0.05))

  # 3: Check that input ind_deg is a list of at least two data.frames.
  # stop when not.
  expect_error(metaRNAseq(ind_deg$DEG2_df,
                          test_statistic = "fishercomb",
                          BHth = 0.05))

  # 4: Check that input ind_deg is a list of at least two data.frames.
  # continue when it is.
  expect_no_error(metaRNAseq(ind_deg,
                             test_statistic = "fishercomb",
                             BHth = 0.05))
})

#commented to reduce check time
#test_that("Annotate_SPID", {
#  websiteLive <- getOption("enrichR.live", default = FALSE)
#  # 1: Check for connectivity
#  if (websiteLive) {
#    expect_no_error(Annotate_SPID(ind_deg$DEG2_df, "WikiPathway_2023_Human"))
#  } else {
#    expect_error(Annotate_SPID(ind_deg$DEG2_df, "WikiPathway_2023_Human"))
#  }
#})

#commented to reduce check time
#test_that("combine_fisher_invnorm", {
#  invnorm <- metaRNAseq(ind_deg, test_statistic = "invnorm", BHth = 0.05, nrep = c(2,2))
#  fishercomb <- metaRNAseq(ind_deg, test_statistic = "fishercomb", BHth = 0.05)
#  # 1: Check the output on toy data
#  expect_no_error(combine_fisher_invnorm(ind_deg,
#                                         invnorm, fishercomb,
#                                         adjpval = 0.05,
#                                         output_tsv = FALSE))
#
#})

#commented to reduce check time
#test_that("DownloadArchS4", {
#  # 1: Check that downloaded counts matrix matches the query
#  GSM <- "GSM3447008"
#  GEO_count_matrix <- DownloadArchS4(GSM, species = "human",
#                                     print_tsv = FALSE, filename = NULL)
#  expect_equal(colnames(GEO_count_matrix), "GSM3447008")
#})

#commented to reduce check time
#test_that("GEOmetadata", {
#  # 1: Check that downloaded metadata matches the query
#  GSE <- "GSE121810"
#  meta <- GEOmetadata(GSE)
#  expect_equal(unique(meta$GSE), "GSE121810")
#})

test_that("SVenn", {
  # 1: Expect no error when the number of sets is up to 7 sets
  S_list <- list(SP1 <- c("EPCAM", "CD24",  "DLK1",  "CDCP1", "LYVE1","TOP2A"),
                 SP2 <- c("DLK1", "EPCAM", "EGFR", "UPK1A", "UPK2"))
  expect_no_error(SVenn(S_list,
                        output_intersectionFile = FALSE,
                        cols.use = c("pink", "yellow")))

  # 2: Expect error when the number of sets is higher than 7
  S_list <- list(SP1 <- c("EPCAM", "CD24",  "DLK1",  "CDCP1", "LYVE1","TOP2A"),
                 SP2 <- c("DLK1", "EPCAM", "EGFR", "UPK1A", "UPK2"),
                 SP3 <- c("MKI67","CD4","CD8"),
                 SP4 <- c("MKI67","CD4","CD8"),
                 SP5 <- c("SELP"),
                 SP6 <- c("SELE"),
                 SP7 <- c("CECAM1"),
                 SP8 <- c("EPCAM", "CD24"))
  expect_error(SVenn(S_list, output_intersectionFile = FALSE))
})

test_that("Splot", {
  # 1: Check the output for specific surface protein genes
  SurfaceProteins_df <- Gene2SProtein(ind_deg$DEG2_df$GeneID, input_type = "gene_name")
  expect_no_error(Splot(SurfaceProteins_df))
})

#commented to reduce check time
#test_that("TCGA_download", {
#GBM_list_s1 <- TCGA_download(project="TCGA-GBM",
#                             whichcounts = "unstranded",
#                             save.matrix = FALSE, save.metadata = FALSE,
#                             barcodes = c("TCGA-06-0878-01A-01R-1849-01"))
#expect_equal(colnames(GBM_list_s1[[1]]), "TCGA-06-0878-01A-01R-1849-01")
#})

#commented to reduce check time
#test_that("Enrichment", {
#  websiteLive <- getOption("enrichR.live", default = FALSE)
#  # 1: Check for connectivity
#  if (!websiteLive) {
#    message("EnrichR website is not reachable. Skipping test.")
#  } else {
#    result <- tryCatch({
#      Enrichment(ind_deg,
#                 enrich.databases = c("GO_Cellular_Component_2021"),
#                 save.results = FALSE)
#    }, error = function(e) {
#      message(paste("Error in Enrichment:", e$message))
#      return(NULL)
#    })
#    # 2: Check successful execution)
#    expect_true(!is.null(result))
#  }
#})

#commented to reduce check time
#test_that("plotPCA", {
#  expect_no_error(SurfR::plotPCA(matrix = cpm(countData),
#          metadata = metadata,
#          nTOP = 100,
#          dims = c(1,2),
#          color.by = "condition", shape.by = "therapy",
#          label = FALSE, main = "PCA"))
#})

#commented to reduce check time
#test_that("DGE", {
#  expect_no_error(DGE(expression = countData, metadata = metadata,
#                  Nreplica = 2,
#                  design = "~condition",condition = "condition",
#                  TEST = "A", CTRL = "B"))
#})
