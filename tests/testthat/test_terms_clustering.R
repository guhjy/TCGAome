context("Tests terms_clustering.R")
source("helper_functions.R")


################################################################################
## Constructor tests
################################################################################

test_that("TermsClustering", {
    random_enrichment <- get_random_enrichment(1000, 10)
    significance_threshold <- get_significance_threshold(
        random_enrichment@raw_enrichment$pvalue,
        10
    )
    adj_method <- "none"
    for (distance_measure in
         random_enrichment@gene_annotations@func_similarity_methods) {

        random_clustering <- TCGAome::TermsClustering(
            random_enrichment,
            distance_measure = distance_measure,
            significance_threshold = significance_threshold,
            adj_method = adj_method)
        testthat::expect_is(random_clustering, "TermsClustering")
        testthat::expect_equal(
            random_clustering@gene_list_enrichment,
            random_enrichment)
        testthat::expect_equal(
            random_clustering@significance_threshold,
            significance_threshold)
        testthat::expect_equal(
            random_clustering@adj_method,
            adj_method)
        expected_significant_results <- TCGAome::get_significant_results(
            random_enrichment,
            significance_threshold = significance_threshold,
            adj_method = adj_method)
        row.names(expected_significant_results) <-
            seq(length=nrow(expected_significant_results))
        testthat::expect_equal(
            random_clustering@significant_results[
                , c("Term", "pvalue", "adj_pvalue")],
            expected_significant_results
            )
        testthat::expect_equal(
            random_clustering@distance_measure,
            distance_measure)
        testthat::expect_is(
            random_clustering@distance_matrix,
            "dist")
    }
})

test_that("TermsClustering invalid arguments", {
    random_enrichment <- get_random_enrichment(1000, 10)
    significance_threshold <- get_significance_threshold(
        random_enrichment@raw_enrichment$pvalue,
        10
    )
    distance_measure = "binary"
    adj_method = "none"
    testthat::expect_error(TCGAome::TermsClustering(
        NULL,
        distance_measure = distance_measure,
        significance_threshold = significance_threshold,
        adj_method = adj_method))
    testthat::expect_error(TCGAome::TermsClustering(
        random_enrichment,
        distance_measure = "non-supported",
        significance_threshold = significance_threshold,
        adj_method = adj_method))
    testthat::expect_error(TCGAome::TermsClustering(
        random_enrichment,
        distance_measure = distance_measure,
        significance_threshold = 1.1,
        adj_method = adj_method))
    testthat::expect_error(TCGAome::TermsClustering(
        random_enrichment,
        distance_measure = distance_measure,
        significance_threshold = -0.1,
        adj_method = adj_method))
    testthat::expect_error(TCGAome::TermsClustering(
        random_enrichment,
        distance_measure = distance_measure,
        significance_threshold = significance_threshold,
        adj_method = "non-supported"))
})

test_that("TermsClustering not enough significant terms", {
    random_enrichment <- get_random_enrichment(1000, 10)
    significance_threshold <- get_significance_threshold(
        random_enrichment@raw_enrichment$pvalue,
        1
    )
    distance_measure = "binary"
    adj_method = "none"
    # No significant terms
    testthat::expect_error(TCGAome::TermsClustering(
        random_enrichment,
        distance_measure = distance_measure,
        significance_threshold = significance_threshold / 2,
        adj_method = adj_method))
    # Only one significant term
    random_enrichment@raw_enrichment[1, "pvalue"] <- significance_threshold / 2
    testthat::expect_error(TCGAome::TermsClustering(
        random_enrichment,
        distance_measure = distance_measure,
        significance_threshold = significance_threshold / 2,
        adj_method = adj_method))
})



