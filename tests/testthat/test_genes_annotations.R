context("Tests genes_annotations.R functions")

## Retrieves a random term from the annotations
get_random_term <- function(gene_annotations) {
    random_term <- gene_annotations@raw_annotations$Term[runif(
        1,
        max=length(gene_annotations@raw_annotations$Term))]
    return(random_term)
}

## Returns a random term with the given frequency
get_term_by_freq <- function(gene_annotations, freq) {
    term_by_freq <- gene_annotations@term2gene[length(unlist(gene_annotations@term2gene$Gene)), ]$Term
    random_term <- term_by_freq[runif(1,max = length(term_by_freq))]
    return(random_term)
}

## Test all term-level methods
test_annotations_methods <- function(gene_annotations) {

    ## Gets term frequency
    random_term <- get_random_term(gene_annotations)
    freq <- TCGAome::get_term_freq(
        gene_annotations,
        random_term)
    testthat::expect_is(freq, "numeric")
    #random_term2 <- get_term_by_freq(gene_annotations, freq)
    #freq2 <- TCGAome::get_term_freq(
    #    gene_annotations,
    #    random_term2)
    #testthat::expect_is(freq2, "numeric")
    #testthat::expect_equal(freq, freq2)

    ## Gets UI similarity
    random_term1 <- get_random_term(gene_annotations)
    random_term2 <- get_random_term(gene_annotations)
    similarity = TCGAome::get_functional_similarity(
        gene_annotations,
        term1 = random_term1,
        term2 = random_term1,
        method = "UI")
    testthat::expect_is(similarity, "numeric")

    ## Gets binary similarity
    random_term1 <- get_random_term(gene_annotations)
    random_term2 <- get_random_term(gene_annotations)
    similarity = TCGAome::get_functional_similarity(
        gene_annotations,
        term1 = random_term1,
        term2 = random_term1,
        method = "binary")
    testthat::expect_is(similarity, "numeric")

    ## Gets cosine similarity
    random_term1 <- get_random_term(gene_annotations)
    random_term2 <- get_random_term(gene_annotations)
    similarity = TCGAome::get_functional_similarity(
        gene_annotations,
        term1 = random_term1,
        term2 = random_term1,
        method = "cosine")
    testthat::expect_is(similarity, "numeric")

    ## Gets Bray-Curtis similarity
    random_term1 <- get_random_term(gene_annotations)
    random_term2 <- get_random_term(gene_annotations)
    similarity = TCGAome::get_functional_similarity(
        gene_annotations,
        term1 = random_term1,
        term2 = random_term1,
        method = "bray-curtis")
    testthat::expect_is(similarity, "numeric")
}

test_annotations_slots <- function(gene_annotations) {
    testthat::expect_is(gene_annotations, "GeneAnnotations")
    testthat::expect_is(gene_annotations@raw_annotations, "data.frame")
    testthat::expect_equal(colnames(gene_annotations@raw_annotations), c("Gene", "Term"))
    testthat::expect_is(gene_annotations@gene2term, "data.frame")
    testthat::expect_is(gene_annotations@term2gene, "data.frame")
    testthat::expect_is(gene_annotations@max_term_freq, "integer")
}

## Test GOA-human annotations
test_that("load_goa()", {
    ## Loads human GOA for a first time
    human_goa = TCGAome::load_goa()
    testthat::expect_is(human_goa, "GeneAnnotations")
    testthat::expect_equal(colnames(human_goa@raw_annotations), c("Gene", "Term"))
    ## Loads the cached version
    human_goa = TCGAome::load_goa()
    testthat::expect_is(human_goa, "GeneAnnotations")
    testthat::expect_equal(colnames(human_goa@raw_annotations), c("Gene", "Term"))
    ## Test methods
    test_annotations_methods(human_goa)
    ## Test slots
    test_annotations_slots(human_goa)
})

## Test GOA-Uniprot annotations
test_that("load_goa(uniprot)", {
    ## Loads human GOA for a first time
    uniprot_goa = TCGAome::load_goa(search_universe = "uniprot")
    testthat::expect_is(uniprot_goa, "GeneAnnotations")
    testthat::expect_equal(colnames(uniprot_goa@raw_annotations), c("Gene", "Term"))
    ## Loads the cached version
    uniprot_goa = TCGAome::load_goa(search_universe = "uniprot")
    testthat::expect_is(uniprot_goa, "GeneAnnotations")
    testthat::expect_equal(colnames(uniprot_goa@raw_annotations), c("Gene", "Term"))
    ## Test methods
    test_annotations_methods(uniprot_goa)
    ## Test slots
    test_annotations_slots(uniprot_goa)
})

test_that("load_goa(other)", {
    ## Loads a GOA for a non supported universe
    other_goa = TCGAome::load_goa(search_universe = "other")
    testthat::expect_null(other_goa)
})

## Test HPO annotations
test_that("load_hpo()", {
    ## Loads HPO for a first time
    hpo = TCGAome::load_hpo()
    testthat::expect_is(hpo, "GeneAnnotations")
    testthat::expect_equal(colnames(hpo@raw_annotations), c("Gene", "Term"))
    ## Loads the cached version
    hpo = TCGAome::load_hpo()
    testthat::expect_is(hpo, "GeneAnnotations")
    testthat::expect_equal(colnames(hpo@raw_annotations), c("Gene", "Term"))
    ## Test methods
    test_annotations_methods(hpo)
    ## Test slots
    test_annotations_slots(hpo)
})

## Test KEGG annotations
test_that("load_kegg()", {
    ## Loads KEGG for a first time
    kegg = TCGAome::load_kegg()
    testthat::expect_is(kegg, "GeneAnnotations")
    testthat::expect_equal(colnames(kegg@raw_annotations), c("Gene", "Term"))
    ## Loads the cached version
    kegg = TCGAome::load_kegg()
    testthat::expect_is(kegg, "GeneAnnotations")
    testthat::expect_equal(colnames(kegg@raw_annotations), c("Gene", "Term"))
    ## Test methods
    test_annotations_methods(kegg)
    ## Test slots
    test_annotations_slots(kegg)
})

## Test OMIM annotations
test_that("load_omim()", {
    ## Loads omim for a first time
    omim = TCGAome::load_omim()
    testthat::expect_is(omim, "GeneAnnotations")
    testthat::expect_equal(colnames(omim@raw_annotations), c("Gene", "Term"))
    ## Loads the cached version
    omim = TCGAome::load_omim()
    testthat::expect_is(omim, "GeneAnnotations")
    testthat::expect_equal(colnames(omim@raw_annotations), c("Gene", "Term"))
    ## Test methods
    test_annotations_methods(omim)
    ## Test slots
    test_annotations_slots(omim)
})
