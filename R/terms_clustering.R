#' @include gene_list_enrichment.R
NULL

setClass("TermsClustering",
         representation(gene_list_enrichment = "GeneListEnrichment",
                        significance_thr = "numeric",
                        adj_method = "character",
                        significant_results = "data.frame",
                        distance_measure = "character",
                        distance_matrix = "dist"),
         prototype(
             gene_list_enrichment = NULL,
             distance_measure = NULL,
             significance_thr = 0.05,
             adj_method = "none"),
         validity = function(object) {
             stopifnot(! is.null(gene_list_enrichment) &&
                           ! is.null(distance_measure) &&
                           ! is.null(significance_thr) &&
                           ! is.null(adj_method))
         }
)

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
                   significance_thr, adj_method){

              stopifnot(distance_measure %in% c("UI", "binary", "bray-curtis", "cosine"))

              ## Initialize input data
              .Object@gene_list_enrichment <- gene_list_enrichment
              .Object@distance_measure <- distance_measure

              ## Gets only significant enrichment results
              significant_results <- get_significant_results(
                  gene_list_enrichment, significance_thr, adj_method)

              ## Computes distance matrix
              .Object@distance_matrix <- get_term_distance_matrix(
                  gene_list_enrichment@gene_annotations, distance_measure,
                  subset = significant_results$Term)

              ## Computes clustering
              clusters <- TCGAome::.pam_clustering(
                  .Object@distance_matrix)
              significant_results <- merge(significant_results, clusters, by = "Term")

              ## Computes MDS
              mds <- TCGAome::.multidimensional_scaling(
                  .Object@distance_matrix)
              significant_results <- merge(significant_results, mds, by = "Term")

              ## Retrieves the frequency of annotation
              significant_results$Freq <- sapply(
                  significant_results$Term,
                  FUN = function(term) {
                      TCGAome::get_term_freq(
                      gene_list_enrichment@gene_annotations, term)
                      }
                  )

              .Object@significant_results <- significant_results

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
