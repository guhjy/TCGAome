context("Tests biomart.R functions")

test_that("load_biomart()", {
  biomart1 = TCGAome::load_biomart()
  expect_is(biomart1, "Mart")
  biomart2 = TCGAome::load_biomart()
  expect_is(biomart2, "Mart")
  expect_equal(biomart1, biomart2)
})


test_that("query_biomart()", {
  results1 = TCGAome::query_biomart(attributes = c("entrezgene"), filters = c("hgnc_symbol"), values = c("BRCA1"))
  expect_equal(dim(results1), c(1,1))
  results2 = TCGAome::query_biomart(attributes = c("entrezgene"), filters = c("hgnc_symbol"), values = c("BRCA1", "CFTR"))
  expect_equal(dim(results2), c(2,1))
})
