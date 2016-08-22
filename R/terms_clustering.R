#' @include gene_list_enrichment.R
NULL


#' TermsClustering class
#'
#' An S4 class to represent the dimensionality reduction computed on
#' enrichment results. Given an object of class "GeneListEnrichment" and a
#' distance measure this class computes:
#' (1) the Partitioning Around Medoids (PAM) clustering
#' (2) the cluster representative selection based on the enrichment significance
#' and the frequency of annotation of the term to give priority to the most
#' specific terms
#' (3) the MultiDimensional Scaling (MDS) results on the cluster representatives
#'
#' @slot gene_list_enrichment The object of class GeneListEnrichment that
#' that contains the enrichment results
#' @slot significance_threshold The significance threshold applied to the enrichment
#' results. The significance threshold is a fixed value, if you want to modify
#' it you would need to create another object.
#' @slot adj_method The multiple test adjustment method. Methods supported are
#' those supported by p.adjust()
#' @slot significant_results A data.frame with the selected annotations terms
#' after applying the multiple test correction and the significance threshold.
#' @slot distance_measure The distance measure employed for the clustering
#' @slot distance_matrix The distance matrix between annotation terms by using
#' the given distance_measure
#'
#' @export
#'
setClass("TermsClustering",
         representation(gene_list_enrichment = "GeneListEnrichment",
                        significance_threshold = "numeric",
                        adj_method = "character",
                        significant_results = "data.frame",
                        distance_measure = "character",
                        distance_matrix = "dist"),
         prototype(
             gene_list_enrichment = NULL,
             distance_measure = NULL,
             significance_threshold = 0.05,
             adj_method = "none"),
         validity = function(object) {
             errors <- character()

             # Check that there are at least 2 terms
             func_similarity_methods <-
                 object@gene_list_enrichment@gene_annotations@func_similarity_methods
             if (! object@distance_measure %in% func_similarity_methods) {
                 msg <- paste(
                     "Not supported distance_measure",
                     object@distance_measure, ". Use one of ",
                     str(func_similarity_methods),
                     sep="")
                 errors <- c(errors, msg)
             }
             if (is.null(object@significance_threshold) ||
                 object@significance_threshold < 0 ||
                 object@significance_threshold > 1) {
                 msg <- paste("Not supported significance_threshold", object@significance_threshold, sep="")
                 errors <- c(errors, msg)
             }
             if (! object@adj_method %in% p.adjust.methods) {
                 msg <- paste("Not supported adj_method ", object@adj_method,
                              ". Use one of ", str(p.adjust.methods), sep="")
                 errors <- c(errors, msg)
             }

             if (length(errors) == 0) TRUE else errors
         }
)

#' TermsClustering constructor
#'
#' The constructor for TermsClustering
#'
#' @param gene_list_enrichment The object of class GeneListEnrichment that
#' that contains the enrichment results
#' @param distance_measure The distance measure employed for the clustering
#' @param significance_threshold The significance threshold applied to the enrichment
#' results. The significance threshold is a fixed value, if you want to modify
#' it you would need to create another object.
#' @param adj_method The multiple test adjustment method. Methods supported are
#' those supported by p.adjust()
#'
#' @return The TermsClustering object
#'
#' @export
#'
#' @examples
#' TODO
TermsClustering <- function(...) new("TermsClustering",...)

# Hierarchical clustering
#hclust_clustering <- function(distance, distance_threshold) {
#    # Clustering with hclust and a static distance threshold
#    clustering <- hclust(distance)
#    clusters <- cutree(clustering, h = distance_threshold)
#    clusters
#}

## PAM and silhouette clustering
.pam_clustering <- function(distance) {
    max_clusters <- attr(distance, "Size") - 1
    ## Finds the optimal number of clusters by silhouette analysis
    sil_width <- sapply(2:max_clusters, FUN = function(x) cluster::pam(distance, x)$silinfo$avg.width)
    names(sil_width) <- 2:max_clusters
    k_best <- as.numeric(names(which.max(sil_width)))
    message(paste("Optimal number of clusters: ", k_best))
    ## Clusters with PAM
    clustering <- cluster::pam(distance, k_best)
    return(data.frame(Term = names(clustering$clustering),
                      Cluster = as.vector(clustering$clustering),
                      stringsAsFactors = F))
}

## Multidimensional Scaling
.multidimensional_scaling <- function(distance) {
    number_elements <- attr(distance, "Size")
    if (number_elements < 2) {
        x <- rep(0, number_elements)
        y <- rep(0, number_elements)
    } else if (number_elements == 2) {
        x <- c(-distance[1] / 2, distance[1] / 2)
        y <- rep(0, number_elements)
    } else {
        fit <- cmdscale(distance, eig = TRUE, k = 2)
        x <- fit$points[, 1]
        if (dim(fit$points)[2] > 1) {
            y <- fit$points[, 2]
        } else {
            y <- rep(0, length(x))
            names(y) <- names(x)
        }
    }
    return(data.frame(Term = names(x), x = x, y = y,
                      stringsAsFactors = F))
}

setMethod("initialize",
          signature(.Object = "TermsClustering"),
          function(.Object, gene_list_enrichment, distance_measure,
                   significance_threshold, adj_method){

              ## Initialize input data
              .Object@gene_list_enrichment <- gene_list_enrichment
              .Object@distance_measure <- distance_measure
              .Object@significance_threshold <- significance_threshold
              .Object@adj_method <- adj_method

              ## Gets only significant enrichment results
              .Object@significant_results <- TCGAome::get_significant_results(
                  gene_list_enrichment, significance_threshold, adj_method)

              ## Computes distance matrix
              .Object@distance_matrix <- TCGAome::get_term_distance_matrix(
                  gene_list_enrichment@gene_annotations, distance_measure,
                  terms_subset = .Object@significant_results$Term)

              ## Checks object validity
              validObject(.Object)

              ## Computes clustering
              clusters <- TCGAome::.pam_clustering(
                  .Object@distance_matrix)
              .Object@significant_results <- merge(
                  .Object@significant_results,
                  clusters,
                  by = "Term")

              ## Computes MDS
              mds <- TCGAome::.multidimensional_scaling(
                  .Object@distance_matrix)
              .Object@significant_results <- merge(
                  .Object@significant_results,
                  mds,
                  by = "Term")

              ## Retrieves the frequency of annotation
              .Object@significant_results$Freq <- sapply(
                  .Object@significant_results$Term,
                  FUN = function(term) {
                      TCGAome::get_term_freq(
                      gene_list_enrichment@gene_annotations, term)
                      }
                  )

              return(.Object)
          })

setGeneric("plot_scatter",
           signature = c("x"),
           function(x) standardGeneric("plot_scatter"))

setMethod("plot_scatter",
          c("x" = "TermsClustering"),
          function(x) {
              NULL
          })
