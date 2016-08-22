## Retrieves a random term from the annotations
get_random_term <- function(gene_annotations) {
    random_term <- gene_annotations@raw_annotations$Term[runif(
        1,
        max=length(gene_annotations@raw_annotations$Term))]
    return(random_term)
}

## Returns a random term with the given frequency
get_term_by_freq <- function(gene_annotations, freq) {
    term_by_freq <- gene_annotations@term2gene[
        sapply(
            gene_annotations@term2gene$Gene,
            FUN= function(x)
                length(x) / gene_annotations@max_term_annotations) == freq,
        ]$Term
    random_term <- term_by_freq[runif(1, min = 1, max = length(term_by_freq))]
    return(random_term)
}

# Returns a random data.frame of annotations
get_random_annotations <- function(size) {
    num_genes <- as.integer(size/10)
    genes <- sapply(1:num_genes,
                    FUN=function(x) paste(
                        sample(c(LETTERS, 0:9),
                               size = 4,
                               replace = TRUE),
                        collapse = ""))
    num_terms <- as.integer(size/10)
    terms <- sapply(1:num_terms,
                    FUN=function(x) paste(
                        "GO:",
                        paste(sample(c(0:9),
                                     size = 4,
                                     replace = TRUE),
                              collapse = ""),
                        sep=""))
    return(data.frame(
        Gene = sample(genes, size = size, replace = TRUE),
        Term = sample(terms, size = size, replace = TRUE)))
}

# Returns a random GeneListEnrichment object
get_random_enrichment <- function(size, size_random_gene_list) {
    random_annotations <- get_random_annotations(size)
    random <- TCGAome::GeneAnnotations(random_annotations, "random")
    random_genes <- sample(random@gene2term$Gene, 10)
    random_enrichment <- TCGAome::GeneListEnrichment(
        gene_list = random_genes,
        gene_annotations = random)
}

# Returns the significance threshold for n most significant terms
get_significance_threshold <- function(pvalues, n) {
    max(pvalues[order(pvalues)[1:n]])
}
