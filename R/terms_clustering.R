#' @include gene_list_enrichment.R
NULL

setClass("TermsClustering",
         representation(gene_list_enrichment = "GeneListEnrichment",
                        significance_thr = "numeric",
                        adj_method = "character",
                        significant_results = "data.frame",
                        distance_measure = "character",
                        distance_matrix = "dist",
                        clusters = "data.frame",
                        mds = "data.frame"),
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

# PAM and silhouette clustering
.pam_clustering <- function(distance) {
    # Clustering with PAM and silhouette
    max_clusters <- attr(distance, "Size") - 1
    sil_width <- numeric(max_clusters)
    #TODO: vectorize this!
    for (k in 2:max_clusters) {
        sil_width[k] <- cluster::pam(distance, k)$silinfo$avg.width
    }
    k_best <- which.max(sil_width[2:max_clusters])
    message(paste("Optimal number of clusters: ", k_best))
    clustering <- cluster::pam(distance, k_best)
    return(data.frame(Term = names(clustering$clustering),
                      Cluster = as.vector(clustering$clustering)))
}

# MDS
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
    return(data.frame(x = x, y = y))
}

setMethod("initialize",
          signature(.Object = "TermsClustering"),
          function(.Object, gene_list_enrichment, distance_measure,
                   significance_thr, adj_method){

              stopifnot(distance_measure %in% c("UI", "binary", "bray-curtis", "cosine"))

              ## Initialize input data
              .Object@gene_list_enrichment <- gene_list_enrichment
              .Object@distance_measure <- distance_measure
              .Object@significant_results <- get_significant_results(
                  gene_list_enrichment, significance_thr, adj_method)
              .Object@distance_matrix <- get_term_distance_matrix(
                  gene_list_enrichment@gene_annotations, distance_measure,
                  subset=.Object@significant_results$Term)

              ## Computes clustering
              .Object@clusters <- TCGAome::.pam_clustering(
                  .Object@distance_matrix)

              ## Computes MDS
              .Object@mds <- TCGAome::.multidimensional_scaling(
                  .Object@distance_matrix)

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
