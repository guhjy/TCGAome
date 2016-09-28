#' @include gene_annotations.R
NULL

#' GeneListEnrichment class
#'
#' An S4 class to represent gene list enrichment results. This class depends on
#' a specific set of gene annotations represented by the class "GeneAnnotations".
#' The enrichment calculated is based on the Fisher exact test and we give
#' support to the usual multiple test adjustment methods.
#'
#' @slot gene_list The list of genes over which the enrichment is computed.
#' @slot raw_enrichment A data.frame with the enrichment results
#'
#' @export
#' @import fastmatch
#'
setClass("GeneListEnrichment",
         representation(gene_list = "character",
                        raw_enrichment = "data.frame"),
         prototype(
             gene_list = c(),
             raw_enrichment = data.frame()),
         validity = function(object) {
             errors <- character()

             # Check that there are at least 2 terms
             if (length(object@gene_list[!is.na(object@gene_list)]) < 1) {
                 msg <- "Empty gene_list provided"
                 errors <- c(errors, msg)
             }

             if (length(errors) == 0) TRUE else errors
         }
)

#' GeneListEnrichment constructor
#'
#' The constructor for GeneListEnrichment
#'
#' @param gene_list The list of genes over which the enrichment is computed.
#' @param gene_annotations The object GeneAnnotations to which the enrichment
#' refers.
#'
#' @return The GeneAnnotations object
#'
#' @export
#'
#' @examples
#' TODO
GeneListEnrichment <- function(...) new("GeneListEnrichment",...)

## Computes the enrichment with the Fisher test for every term
.compute_enrichment <- function(gene_list, gene_annotations) {
    require(fastmatch)
    all_genes <- gene_annotations@gene2term$Gene
    all_terms <- gene_annotations@term2gene$Term
    pvalues <- as.double(
        vapply(all_terms,
               FUN = function(term) {
                   associated_genes = as.character(unlist(
                       gene_annotations@term2gene[term, ]$Gene))
                   intersection_genes <- gene_list[!is.na(fmatch(gene_list, associated_genes))]
                   union_genes <- c(gene_list, associated_genes[is.na(fmatch(associated_genes, gene_list))])
                   fisher.test(matrix(c(
                       # intersection between gene_list and associated_genes
                       length(intersection_genes),
                       # associated_genes not in gene_list
                       length(associated_genes [is.na(fmatch(associated_genes, intersection_genes))]),
                       # gene_list not in associated_genes
                       length(gene_list [is.na(fmatch(gene_list, intersection_genes))]),
                       # genes not in associated_genes and not in gene_list
                       length(all_genes [is.na(fmatch(all_genes, union_genes))])
                    ), nrow = 2, ncol = 2),
                    alternative = "greater")$p.value
                },
               FUN.VALUE = numeric(1)))
    return(pvalues)
}

setMethod("initialize",
          signature(.Object = "GeneListEnrichment"),
          function(.Object, gene_list, gene_annotations){
              ## Initialize input data
              futile.logger::flog.info("Initializing GeneListEnrichment...")
              .Object@gene_list <- gene_list
              futile.logger::flog.info("Input gene list is of length : %d", length(.Object@gene_list))

              ## Checks object validity
              futile.logger::flog.info("Checking object's validity...")
              validObject(.Object)
              futile.logger::flog.info("Object is valid!")

              ## Clear NA_character genes
              .Object@gene_list <- .Object@gene_list[!is.na(.Object@gene_list)]
              futile.logger::flog.info("After removing missing values input gene list is of length : %d", length(.Object@gene_list))

              ## Computes enrichment
              futile.logger::flog.info("Computing Fisher test enrichment on the gene list...")
              pvalues <- .compute_enrichment(.Object@gene_list, gene_annotations)
              futile.logger::flog.info("Most significant p-value %s and less significant %s", min(pvalues), max(pvalues))

              ## Computes the relative frequency of every term
              #freqs <- as.double(sapply(all_terms, function(term){
              #    TCGAome::get_term_freq(gene_annotations, term)
              #}))
              ## Builds results data.frame
              all_terms <- gene_annotations@term2gene$Term
              .Object@raw_enrichment <- data.frame(
                  Term = all_terms,
                  pvalue = pvalues,
                  row.names = NULL, stringsAsFactors = F)

              return(.Object)
          })


#' get_significant_results()
#'
#' Function that gets the significant results from the enrichment results.
#'
#' @param x The GeneListEnrichment class on which the method will run.
#' @param significance_thr The significance enrichment threshold
#' @param adj_method The multiple test correction method (any of those supported
#' by p.adjust())

#' @return A data.frame with the significant results
#'
#' @export
#'
#' @examples
#' TODO
setGeneric("get_significant_results",
           signature = c("x", "significance_threshold", "adj_method"),
           function(x, significance_threshold, adj_method)
               standardGeneric("get_significant_results"))

#' @aliases get_significant_results
#' @export
setMethod("get_significant_results",
          c("x" = "GeneListEnrichment", "significance_threshold" = "numeric",
            "adj_method" = "character"),
          function(x, significance_threshold, adj_method) {
              if (significance_threshold < 0 || significance_threshold > 1) {
                  stop(paste("Not valid significance threshold '", significance_threshold, "'. It must be in the range [0, 1]", sep = ""))
              }
              if (! adj_method %in% p.adjust.methods) {
                  stop(paste("Non supported p-value adjustment method: '", adj_method, "'", sep=""))
              }
              x@raw_enrichment$adj_pvalue <- p.adjust(x@raw_enrichment$pvalue, method = adj_method)
              return(x@raw_enrichment[x@raw_enrichment$adj_pvalue <= significance_threshold, ])
          })

#' get_terms_clustering()
#'
#' Function that creates a new object of class TermsClustering. The clustering
#' is performed on a significant set of enrichment results that is defined by
#' the parameters passed to this function.
#'
#' @param x The GeneListEnrichment class on which the method will run.
#' @param gene_annotations The GeneAnnotations on which the enrichment was computed.
#' @param distance_measure The similarity distance metric employed.
#' @param significance_threshold The significance enrichment threshold
#' @param adj_method The multiple test correction method (any of those supported
#' by p.adjust())

#' @return An object of class TermsClustering
#'
#' @export
#'
#' @examples
#' TODO
setGeneric("get_terms_clustering",
           signature = c("x", "gene_annotations", "distance_measure",
                         "significance_threshold",
                         "adj_method",
                         "max_clusters"),
           function(x, gene_annotations, distance_measure,
                    significance_threshold, adj_method, max_clusters)
               standardGeneric("get_terms_clustering"))

#' @aliases get_terms_clustering
#' @export
setMethod("get_terms_clustering", c("x" = "GeneListEnrichment",
                                    "gene_annotations" = "GeneAnnotations",
                                    "distance_measure" = "character",
                                    "significance_threshold" = "numeric",
                                    "adj_method" = "character",
                                    "max_clusters" = "numeric"),
          function(x, gene_annotations, distance_measure,
                   significance_threshold, adj_method, max_clusters) {
              return(TCGAome::TermsClustering(
                  gene_annotations = gene_annotations,
                  gene_list_enrichment = x,
                  distance_measure = distance_measure,
                  significance_thr = significance_threshold,
                  adj_method = adj_method,
                  max_clusters = max_clusters)
              )
          })

