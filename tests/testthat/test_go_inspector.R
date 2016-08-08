context("Tests go_inspector.R functions")

test_that("load_go_ontology()", {
    go_tree1 = TCGAome::load_go_ontology()
    expect_is(go_tree1, "Ontology")
    go_tree2 = TCGAome::load_go_ontology()
    expect_is(go_tree2, "Ontology")
    expect_equal(go_tree1, go_tree2)
})


#test_that("query_biomart()", {
#    results1 = TCGAome::query_biomart(attributes = c("entrezgene"), filters = c("hgnc_symbol"), values = c("BRCA1"))
#    expect_equal(dim(results1), c(1,1))
#    results2 = TCGAome::query_biomart(attributes = c("entrezgene"), filters = c("hgnc_symbol"), values = c("BRCA1", "CFTR"))
#    expect_equal(dim(results2), c(2,1))
#})
