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
#' @slot significance_threshold The significance threshold applied to the enrichment
#' results. The significance threshold is a fixed value, if you want to modify
#' it you would need to create another object.
#' @slot adj_method The multiple test adjustment method. Methods supported are
#' those supported by p.adjust()
#' @slot max_clusters The maximum number of clusters to evaluate (default: 10)
#' @slot significant_results A data.frame with the selected annotations terms
#' after applying the multiple test correction and the significance threshold.
#' @slot distance_measure The distance measure employed for the clustering
#' @slot distance_matrix The distance matrix between annotation terms by using
#' the given distance_measure
#' @slot explained_variance The percentage of variance explained by each MDS
#' component
#' @slot explained_variance_repr The percentage of variance explained by each
#' MDS component computed only to cluster representatives
#' @slot silhouette The silhouette average width for each number of clusters
#'
#' @export
#'
setClass("TermsClustering",
         representation(distance_measure = "character",
                        significance_threshold = "numeric",
                        adj_method = "character",
                        max_clusters = "numeric",
                        significant_results = "data.frame",
                        distance_matrix = "dist",
                        explained_variance = "data.frame",
                        explained_variance_repr = "data.frame",
                        silhouette = "data.frame"),
         prototype(
             distance_measure = NULL,
             significance_threshold = 0.05,
             adj_method = "none",
             max_clusters = 10),
         validity = function(object) {
             errors <- character()

             # Check that there are at least 2 terms
             #FIXME: fix this
             #func_similarity_methods <-
             #    object@gene_list_enrichment@gene_annotations@func_similarity_methods
             #if (! object@distance_measure %in% func_similarity_methods) {
             #    msg <- paste(
             #        "Not supported distance_measure",
             #        object@distance_measure, ". Use one of ",
             #        str(func_similarity_methods),
             #        sep="")
             #    errors <- c(errors, msg)
             #}
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
             if (! is.null(object@max_clusters) & object@max_clusters < 0) {
                 msg <- "Maximum number of clusters has to be a positive integer"
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
#' @slot max_clusters The maximum number of clusters to evaluate (default: 10)
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
.pam_clustering <- function(distance, max_clusters = 10) {
    ## if the max_clusters is set to null or 0 then it
    ## takes the number of elements minus 1
    if (is.null(max_clusters) | max_clusters == 0) {
        max_clusters <- attr(distance, "Size") - 1
    }
    max_clusters <- min(max_clusters,
                       attr(distance, "Size") - 1)
    ## Finds the optimal number of clusters by silhouette analysis
    sil_width <- sapply(2:max_clusters, FUN = function(x) cluster::pam(distance, x)$silinfo$avg.width)
    names(sil_width) <- 2:max_clusters
    k_best <- as.numeric(names(which.max(sil_width)))
    message(paste("Optimal number of clusters: ", k_best))
    ## Clusters with PAM
    clustering <- cluster::pam(distance, k_best)
    return(list(
        clustering = data.frame(Term = names(clustering$clustering),
                                Cluster = as.vector(clustering$clustering),
                                stringsAsFactors = F),
        silhouette = data.frame(k = as.numeric(names(sil_width)),
                                avg_width = as.numeric(sil_width),
                                stringsAsFactors = F)
            )
    )
}

## Multidimensional Scaling
.multidimensional_scaling <- function(distance) {
    # TODO: return explained variance for each component
    number_elements <- attr(distance, "Size")
    if (number_elements < 2) {
        pc1 <- rep(0, number_elements)
        pc2 <- rep(0, number_elements)
        pc3 <- rep(0, number_elements)
        explained_variance = c(1, 0, 0)
    } else if (number_elements == 2) {
        pc1 <- c(-distance[1] / 2, distance[1] / 2)
        pc2 <- rep(0, number_elements)
        pc3 <- rep(0, number_elements)
        explained_variance = c(1, 0, 0)
    } else {
        fit <- cmdscale(distance, eig = TRUE, k = 3)
        pc1 <- fit$points[, 1]
        if (dim(fit$points)[2] > 1) {
            pc2 <- fit$points[, 2]
        } else {
            pc2 <- rep(0, length(x))
            names(y) <- names(x)
        }
        if (dim(fit$points)[2] > 2) {
            pc3 <- fit$points[, 3]
        } else {
            pc3 <- rep(0, length(x))
            names(z) <- names(x)
        }
        explained_variance = fit$eig / sum(fit$eig)
    }
    return(list(
        mds = data.frame(Term = names(pc1), pc1 = pc1, pc2 = pc2, pc3 = pc3,
                         stringsAsFactors = F),
        explained_variance = data.frame(
            component = 1:length(explained_variance),
            explained_variance = explained_variance,
            stringsAsFactors = F)))
}

## Selects the cluster representatives based on significance and frequency
.select_representatives <- function(significant_results, freq_threshold = 0.05) {
    representatives <- as.data.frame(t(sapply(
        unique(significant_results$Cluster),
        FUN = function(cluster) {
            candidates <- significant_results[
                significant_results$Cluster == cluster
                & significant_results$Freq <= freq_threshold, ]
            if (nrow(candidates) == 0) {
                candidates <- significant_results[
                    significant_results$Cluster == cluster,
                    ]
            }
            list(
                Cluster = cluster,
                repr_term =
                    candidates[
                        which.min(candidates$adj_pvalue), "Term"])
        })))
    return(representatives)
}

setMethod("initialize",
          signature(.Object = "TermsClustering"),
          function(.Object, gene_annotations, gene_list_enrichment, distance_measure,
                   significance_threshold, adj_method, max_clusters){

              ## Initialize input data
              .Object@distance_measure <- distance_measure
              .Object@significance_threshold <- significance_threshold
              .Object@adj_method <- adj_method

              ## Gets only significant enrichment results
              .Object@significant_results <- TCGAome::get_significant_results(
                  gene_list_enrichment, significance_threshold, adj_method)

              ## Computes distance matrix
              .Object@distance_matrix <- TCGAome::get_term_distance_matrix(
                  gene_annotations, distance_measure,
                  terms_subset = .Object@significant_results$Term)

              ## Checks object validity
              validObject(.Object)

              ## Computes clustering
              clusters <- TCGAome::.pam_clustering(
                  .Object@distance_matrix,
                  max_clusters = max_clusters)
              .Object@significant_results <- merge(
                  .Object@significant_results,
                  clusters$clustering,
                  by = "Term")
              .Object@silhouette <- clusters$silhouette

              ## Computes MDS
              mds <- TCGAome::.multidimensional_scaling(
                  .Object@distance_matrix)
              .Object@significant_results <- merge(
                  .Object@significant_results,
                  mds$mds,
                  by = "Term")
              .Object@explained_variance <- mds$explained_variance

              ## Retrieves the frequency of annotation
              .Object@significant_results$Freq <- sapply(
                  .Object@significant_results$Term,
                  FUN = function(term) {
                      TCGAome::get_term_freq(
                          gene_annotations, term)
                  }
              )

              ## Select the representative term for every cluster
              representatives <- .select_representatives(
                  .Object@significant_results)
              .Object@significant_results <- merge(
                  .Object@significant_results,
                  representatives,
                  by = "Cluster")

              #TODO: compute MDS again just on representatives
              distance_matrix_repr <- as.dist(as.matrix(
                  .Object@distance_matrix)[
                      unlist(representatives$repr_term),
                      unlist(representatives$repr_term)])
              mds_repr <- TCGAome::.multidimensional_scaling(
                  distance_matrix_repr)
              names(mds_repr$mds) = c("Term", "pc1_repr",
                                      "pc2_repr", "pc3_repr")

              .Object@significant_results <- merge(
                  .Object@significant_results,
                  mds_repr$mds,
                  by = "Term",
                  all = TRUE)
              .Object@explained_variance_repr <- mds_repr$explained_variance

              return(.Object)
          })

.plot_scatter <- function(significant_results, all, pvalue_threshold) {
    ## Prepares the ellipses per cluster
    if (all) {
        data_ell <- data.frame()
        for(cluster in unique(significant_results$Cluster)){

            data_cluster <- significant_results[
                significant_results$Cluster==cluster,]

            if (length(data_cluster$x) > 1) {
                data_points <- as.data.frame(
                    ellipse::ellipse(
                        cor(data_cluster$x, data_cluster$y),
                        scale=c(sd(data_cluster$x),
                                sd(data_cluster$y)),
                        centre=c(mean(data_cluster$x),
                                 mean(data_cluster$y))))
            } else {
                r <- 0.1
                theta <- sort(2*pi*runif(30))
                data_points <- data.frame(
                    x = c(data_cluster$x + r*sin(theta)),
                    y = c(data_cluster$y + r*cos(theta)))
            }
            data_ell <- rbind(
                data_ell,
                cbind(data_points, Cluster=cluster))
        }
    }

    ## Bins frequency for proper representation
    significant_results$Freq_binned = with(
        significant_results,
        ifelse(Freq < 0.05, 1,
               ifelse(Freq <= 0.1, 2, 3)))

    ## Prepares term datasets
    representatives_data <- significant_results[significant_results$Term == significant_results$repr_term, ]
    others_data <- significant_results[significant_results$Term != significant_results$repr_term, ]

    my_palette <- colorRampPalette(RColorBrewer::brewer.pal(
        3,
        "Spectral"))

    plot <- ggplot2::ggplot(data = representatives_data)

    if (all) {
        plot <- plot +
            # Plot an ellipse around each cluster
            ggplot2::geom_path(
                data = data_ell,
                ggplot2::aes(x = x,
                             y = y,
                             group = Cluster),
                size = 0.2,
                linetype = 2,
                colour = I("blue"))

        # Plots all terms in a cluster other than representatives
        plot <- plot +
            ggplot2::geom_point(
                data = others_data,
                ggplot2::aes(x = x,
                             y = y,
                             colour = adj_pvalue,
                             size = Freq),
            shape = I(21))
    }

    plot <- plot +
        # Plots representative terms
        ggplot2::geom_point(
            data = representatives_data,
            ggplot2::aes(x = x,
                         y = y,
                         colour = adj_pvalue,
                         size = Freq),
            shape = I(16),
            alpha = I(0.85)) +
        # Labels representative terms
        ggplot2::geom_text(
            data = representatives_data,
            ggplot2::aes(x = x,
                         y = y,
                         label = Term),
            colour = I(ggplot2::alpha("black", 0.85)),
            size = 3,
            nudge_y = -0.03,
            vjust = "inward",
            hjust = "inward") +
        # Colours by p-value
        ggplot2::scale_colour_gradientn(
            "Significance",
            colours = my_palette(5),
            limits = c(0, pvalue_threshold)) +
        # Size by frequency of annotation
        ggplot2::scale_size(
            "Frequency",
            range = c(1, 4),
            label = scales::percent,
            #breaks = unique(significant_results$Freq_binned),
            #labels = c("< 5%", "5-10%", "> 10%"),
            #trans = "log",
            guide = ggplot2::guide_legend(override.aes = list(colour = "gray"))) +
        ggplot2::theme_bw()

    plot
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
.multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}


.multiplot_shared_legend <- function(...,
                                     nrow = 1,
                                     ncol = length(list(...)),
                                     position = c("bottom", "right")) {
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplot2::ggplotGrob(
        plots[[1]] + ggplot2::theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) {
        x + ggplot2::theme(legend.position = "none")
    })
    gl <- c(gl, nrow = nrow, ncol = ncol)

    combined <- switch(position,
                       "bottom" = gridExtra::arrangeGrob(
                           do.call(gridExtra::arrangeGrob, gl),
                           legend,
                           ncol = 1,
                           heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = gridExtra::arrangeGrob(
                           do.call(gridExtra::arrangeGrob, gl),
                           legend,
                           ncol = 2,
                           widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth)))
    grid::grid.newpage()
    grid::grid.draw(combined)

}

#' plot_mds()
#'
#' Plots the enriched terms on the three first components of the MDS results.
#'
#' @param x The GeneTermClustering class on which the method will run.
#' @param all A flag indicating if showing all members of each cluster or
#' just the cluster representatives

#' @return The ggplot2 plot
#'
#' @export
#'
#' @examples
#' TCGAome::plot_mds(hpo_term_clustering, all = TRUE)
setGeneric("plot_mds",
           signature = c("x", "all"),
           function(x, all) standardGeneric("plot_mds"))

#' @aliases plot_mds
#' @export
setMethod("plot_mds",
          c("x" = "TermsClustering", "all" = "logical"),
          function(x, all=FALSE) {

              ## Plots components 1 and 2
              if (all) {
                  x@significant_results$y <- x@significant_results$pc1
                  x@significant_results$x <- x@significant_results$pc2
              } else {
                  x@significant_results$y <- x@significant_results$pc1_repr
                  x@significant_results$x <- x@significant_results$pc2_repr
              }
              plot1 = .plot_scatter(x@significant_results, all,
                                    pvalue_threshold = x@significance_threshold)
              plot1 <- plot1 + ggplot2::labs(y = "Comp 1", x = "Comp 2")

              ## Plots components 1 and 3
              if (all) {
                  x@significant_results$y <- x@significant_results$pc1
                  x@significant_results$x <- x@significant_results$pc3
              } else {
                  x@significant_results$y <- x@significant_results$pc1_repr
                  x@significant_results$x <- x@significant_results$pc3_repr
              }
              plot2 = .plot_scatter(x@significant_results, all,
                                    pvalue_threshold = x@significance_threshold)
              plot2 <- plot2 + ggplot2::labs(y = "Comp 1", x = "Comp 3")

              ## Plots components 2 and 3
              if (all) {
                  x@significant_results$x <- x@significant_results$pc2
                  x@significant_results$y <- x@significant_results$pc3
              } else {
                  x@significant_results$x <- x@significant_results$pc2_repr
                  x@significant_results$y <- x@significant_results$pc3_repr
              }
              plot3 = .plot_scatter(x@significant_results, all,
                                    pvalue_threshold = x@significance_threshold)
              plot3 <- plot3 + ggplot2::labs(x = "Comp 2", y = "Comp 3")

              #.multiplot(plot1, plot3, plot2, plot4, cols = 2)
              .multiplot_shared_legend(plot1, plot2, plot3, nrow = 2, ncol = 2, position = "right")
          })

#' plot_explained_variance()
#'
#' Plots the explained variance of each MDS component
#'
#' @param x The GeneTermClustering class on which the method will run.

#' @return The ggplot2 plot
#'
#' @export
#'
#' @examples
#' TCGAome::plot_explained_variance(hpo_term_clustering)
setGeneric("plot_explained_variance",
           signature = c("x"),
           function(x) standardGeneric("plot_explained_variance"))

#' @aliases plot_explained_variance
#' @export
setMethod("plot_explained_variance",
          c("x" = "TermsClustering"),
          function(x) {

              data <- rbind(
                  cbind(set = "all",  x@explained_variance),
                  cbind(set = "representatives",  x@explained_variance_repr))

              ## Plots explained variance
              plot <- ggplot2::ggplot(
                  data = data[data$component <= 10, ],
                  ggplot2::aes(x = component,
                               y = explained_variance,
                               colour = set)) +
                  ggplot2::geom_point(
                      shape = I(16),
                      size = 3) +
                  ggplot2::geom_line(size = 0.2,
                                     linetype = 2) +
                  ggplot2::scale_x_continuous(
                      name = "Component",
                      breaks=data$component) +
                  ggplot2::scale_y_continuous(
                      name = "Explained variance",
                      label = scales::percent) +
                  ggplot2::scale_color_discrete(
                      name = NULL)
                  ggplot2::theme_bw()

              plot
          })

#' plot_silhouette_analysis()
#'
#' Plots the silhouette average width of each number of clusters
#'
#' @param x The GeneTermClustering class on which the method will run.

#' @return The ggplot2 plot
#'
#' @export
#'
#' @examples
#' TCGAome::plot_silhouette_analysis(hpo_term_clustering)
setGeneric("plot_silhouette_analysis",
           signature = c("x"),
           function(x) standardGeneric("plot_silhouette_analysis"))

#' @aliases plot_silhouette_analysis
#' @export
setMethod("plot_silhouette_analysis",
          c("x" = "TermsClustering"),
          function(x) {

              ## Plots silhouette analysis
              plot <- ggplot2::ggplot(
                  data = x@silhouette,
                  ggplot2::aes(x = k, y = avg_width)) +
                  ggplot2::geom_vline(
                      xintercept = x@silhouette[which.max(x@silhouette$avg_width) , ]$k,
                      size = 0.2,
                      linetype = 2,
                      colour = "red") +
                  ggplot2::geom_line(colour="blue",
                                     size = 0.2,
                                     linetype = 2) +
                  ggplot2::geom_point(
                      shape = I(16),
                      size = 3) +
                  ggplot2::scale_x_continuous(
                      name = "Number of clusters",
                      breaks = x@silhouette$k) +
                  ggplot2::scale_y_continuous(
                      name = "Average silhouette width"
                  ) +
                  ggplot2::theme_bw()

              plot
          })
