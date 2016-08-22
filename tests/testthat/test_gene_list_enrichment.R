context("Tests gene_list_enrichment.R")
source("helper_functions.R")


################################################################################
## Constructor tests
################################################################################

test_that("empty gene_list", {
    random_annotations <- get_random_annotations(1000)
    random <- TCGAome::GeneAnnotations(random_annotations, "random")
    empty_genes <- c(NA_character_)
    testthat::expect_error(TCGAome::GeneListEnrichment(
        gene_list = empty_genes,
        gene_annotations = random))
})

test_that("gene_list with NAs", {
    random_annotations <- get_random_annotations(1000)
    random <- TCGAome::GeneAnnotations(random_annotations, "random")
    random_genes_with_NAs <- c(sample(random@gene2term$Gene, 10),
                               c(NA_character_, NA_character_))
    random_enrichment <- TCGAome::GeneListEnrichment(
        gene_list = random_genes_with_NAs,
        gene_annotations = random)
    testthat::expect_is(random_enrichment, "GeneListEnrichment")
    ## the stored gene list does not contain NA values
    testthat::expect_equal(sum(is.na(random_enrichment@gene_list)), 0)
    testthat::expect_equal(sum(is.na(random_genes_with_NAs)), 2)
})

test_that("gene_list not in annotations", {
    random_annotations <- get_random_annotations(1000)
    random <- TCGAome::GeneAnnotations(random_annotations, "random")
    non_annotated_genes <- tolower(sample(random@gene2term$Gene, 10))
    random_enrichment <- TCGAome::GeneListEnrichment(
        gene_list = non_annotated_genes,
        gene_annotations = random)
    testthat::expect_is(random_enrichment, "GeneListEnrichment")
    ## all p-values are equal to 1 as no term can be enriched at all
    testthat::expect_equal(
        random_enrichment@raw_enrichment$pvalue,
        rep(1, length(random_enrichment@raw_enrichment$pvalue)))
})

test_that("random annotations", {
    random_enrichment <- get_random_enrichment(1000, 10)
    random_annotations <- get_random_annotations(1000)
    random <- TCGAome::GeneAnnotations(random_annotations, "random")
    random_genes <- sample(random@gene2term$Gene, 10)
    random_enrichment <- TCGAome::GeneListEnrichment(
        gene_list = random_genes,
        gene_annotations = random)
    testthat::expect_is(random_enrichment, "GeneListEnrichment")
    testthat::expect_equal(random_enrichment@gene_list, random_genes)
    testthat::expect_is(random_enrichment@gene_annotations, "GeneAnnotations")
    testthat::expect_equal(
        colnames(random_enrichment@raw_enrichment),
        c("Term", "pvalue"))
})


################################################################################
## get_significant_results tests
################################################################################

test_that("get_significant_results", {
    random_annotations <- get_random_annotations(1000)
    random <- TCGAome::GeneAnnotations(random_annotations, "random")
    random_genes <- sample(random@gene2term$Gene, 10)
    random_enrichment <- TCGAome::GeneListEnrichment(
        gene_list = random_genes,
        gene_annotations = random)
    significance_threshold_10terms <-
        get_significance_threshold(random_enrichment@raw_enrichment$pvalue,
                                   10)
    significant_results <- TCGAome::get_significant_results(
        random_enrichment,
        significance_threshold = significance_threshold_10terms,
        adj_method = "none")
    testthat::expect_is(significant_results, "data.frame")
    testthat::expect_equal(significant_results$pvalue,
                        significant_results$adj_pvalue)
    testthat::expect_lte(max(significant_results$pvalue),
                         significance_threshold_10terms)
})

test_that("get_significant_results invalid arguments", {
    random_annotations <- get_random_annotations(1000)
    random <- TCGAome::GeneAnnotations(random_annotations, "random")
    random_genes <- sample(random@gene2term$Gene, 10)
    random_enrichment <- TCGAome::GeneListEnrichment(
        gene_list = random_genes,
        gene_annotations = random)
    invalid_significance_threshold <- 1.1
    testthat::expect_error(TCGAome::get_significant_results(
        random_enrichment,
        significance_threshold = invalid_significance_threshold,
        adj_method = "none"))
    invalid_significance_threshold <- -0.1
    testthat::expect_error(TCGAome::get_significant_results(
        random_enrichment,
        significance_threshold = invalid_significance_threshold,
        adj_method = "none"))
    valid_significance_threshold <- 0.1
    invalid_adj_method <- "lacuentalavieja"
    testthat::expect_error(TCGAome::get_significant_results(
        random_enrichment,
        significance_threshold = valid_significance_threshold,
        adj_method = invalid_adj_method))
})

################################################################################
## get_enrichment tests
################################################################################

test_that("get_terms_clustering", {
    random_annotations <- get_random_annotations(1000)
    random <- TCGAome::GeneAnnotations(random_annotations, "random")
    random_genes <- sample(random@gene2term$Gene, 10)
    random_enrichment <- TCGAome::GeneListEnrichment(
        gene_list = random_genes,
        gene_annotations = random)
    significance_threshold_10terms <- get_significance_threshold(
        random_enrichment@raw_enrichment$pvalue,
        10
    )
    terms_clustering <- TCGAome::get_terms_clustering(
        random_enrichment,
        distance_measure = "binary",
        significance_threshold = significance_threshold_10terms,
        adj_method = "none")
    testthat::expect_is(terms_clustering, "TermsClustering")
})

