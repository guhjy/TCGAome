## Loads the Ensembl Biomart for Homo sapiens
load_biomart <- function() {
    ensembl <- attr(load_biomart, "cached_ensembl")
    if (is.null(ensembl)) {
        ## Loads ensembl biomart for Homo sapiens ensembl = useMart('ensembl',dataset='hsapiens_gene_ensembl')
        futile.logger::flog.debug("Loading biomart for the first time")
        ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
        dataset <- "hsapiens_gene_ensembl"
        ensembl <- biomaRt::useDataset(dataset, mart = ensembl_mart)
        attr(load_biomart, "cached_ensembl") <<- ensembl
    } else {
      futile.logger::flog.debug("Returning cached biomart")
    }
    ensembl
}

## Queries the biomart
query_biomart <- function(attributes, filters, values) {
    ## Example mapping HUGO gene names to Entrez gene IDs using BioMart
    biomaRt::getBM(attributes = attributes, filters = filters, values = values, mart = load_biomart())
}
