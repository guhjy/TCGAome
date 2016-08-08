context("Tests pipeline.R functions")

test_that("get_TCGAome_tumor_types()", {
  expect_equal(TCGAome::get_TCGAome_tumor_types(),
               c("ACC", "BLCA", "BRCA", "CESC", "COAD", "COADREAD", "DLBC",
                 "GBM", "HNSC", "KIRC", "KIRP", "LUAD", "LUSC", "OV", "PAAD",
                 "PCPG", "READ", "SKCM", "TGCT", "THYM", "UCEC", "UCS", "UVM"))
})

test_that("get_TCGAome_tumor_types()", {
  expect_equal(TCGAome::get_TCGAome_tumor_types(),
               c("ACC", "BLCA", "BRCA", "CESC", "COAD", "COADREAD", "DLBC",
                 "GBM", "HNSC", "KIRC", "KIRP", "LUAD", "LUSC", "OV", "PAAD",
                 "PCPG", "READ", "SKCM", "TGCT", "THYM", "UCEC", "UCS", "UVM"))
})
