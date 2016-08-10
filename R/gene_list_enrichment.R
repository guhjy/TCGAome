#' @include gene_annotations.R
NULL

setClass("GeneListEnrichment",
         representation(gene_list = "character",
                        gene_annotations = "GeneAnnotations",
                        raw_enrichment = "data.frame"),
         prototype(
             gene_list = c(),
             gene_annotations = NULL),
         validity = function(object) {
             stopifnot(length(object@gene_list[!is.na(object@gene_list)]) > 0 &&
                           ! is.null(object@gene_annotations))
         }
)

setMethod("initialize",
          signature(.Object = "GeneListEnrichment"),
          function(.Object, gene_list, gene_annotations){
              ## Initialize input data
              .Object@gene_list <- gene_list[!is.na(gene_list)]
              .Object@gene_annotations <- gene_annotations
              ## Computes the enrichment with the Fisher test for every term
              all_genes <- gene_annotations@gene2term$Gene
              selected_genes <- .Object@gene_list
              all_terms <- gene_annotations@term2gene$Term
              results <- sapply(all_terms, function(term) {
                  associated_genes = as.character(unlist(
                      gene_annotations@term2gene[gene_annotations@term2gene$Term == term, ]$Gene))
                  fisher.test(matrix(c(
                      length(intersect(selected_genes, associated_genes)),
                      length(associated_genes[! associated_genes %in% selected_genes]),
                      length(selected_genes[! selected_genes %in% associated_genes]),
                      length(all_genes[! all_genes %in% union(associated_genes, selected_genes)])
                  ), nrow = 2, ncol = 2),
                  alternative = "greater")$p.value
              })
              ## Computes the relative frequency of every term
              freqs <- as.double(sapply(all_terms, function(term){
                  TCGAome::get_term_freq(gene_annotations, term)
              }))
              ## Builds results data.frame
              .Object@raw_enrichment <- data.frame(
                  Term = all_terms,
                  pvalue = as.double(results),
                  freq = freqs,
                  row.names = NULL, stringsAsFactors = F)

              return(.Object)
          })

setGeneric("get_significant_results",
           signature = c("x", "significance_thr", "adj_method"),
           function(x, significance_thr, adj_method) standardGeneric("get_significant_results"))

setMethod("get_significant_results",
          c("x" = "GeneListEnrichment", "significance_thr" = "numeric", "adj_method" = "character"),
          function(x, significance_thr, adj_method) {
              stopifnot(significance_thr >= 0 && significance_thr <= 1)
              x@raw_enrichment$adj_pvalue <- p.adjust(x@raw_enrichment$pvalue, method = adj_method)
              return(x@raw_enrichment[x@raw_enrichment$adj_pvalue <= significance_thr, ])
          })
