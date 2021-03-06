% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.r
\name{run_TCGAome}
\alias{run_TCGAome}
\title{Runs the TCGAome analysis pipeline on the given tumor types.}
\usage{
run_TCGAome(tumor_types, run_pca_analysis = TRUE,
  run_hclust_analysis = FALSE, run_rgcca = FALSE, run_rcca = FALSE,
  topN = 5, spls_selection_method = "loadings",
  GO_similarity_measure = "Wang", GOA_search_universe = "human",
  enrichment_significance_thr = 0.01, GO_ontology = "BP",
  multiple_test_adjustment = FALSE)
}
\arguments{
\item{tumor_types}{Vector of tumor types to be analyzed from those available at TCGA (run RTCGAToolbox::getFirehoseDatasets() to see all available types)}

\item{run_pca_analysis}{Flag indicating if the PCA should run, this is part of the data pre-analysis (default: TRUE)}

\item{run_hclust_analysis}{Flag indicating if the hierarchical clustering analysis should run, this is part of the data pre-analysis (default: FALSE)}

\item{run_rgcca}{Flag indicating if the Regularized Generalized Canonical Correlation Analysis (RGCCA) should run (default: FALSE). Beware that this analysis is not intended for datasets with number of variables >> number of samples.}

\item{run_rcca}{Flag indicating if the Regularized Canonical Correlation Analysis (rCCA) should run (default: FALSE). Beware that this analysis is not intended for datasets with number of variables >> number of samples.}

\item{topN}{Indicates the top number of variables to select on MCIA and sPLS results (default: 5). It will select N variables on each of the data types, on each of the three first components and on each extreme of range, that is a maximum of 2*3*2*N, considering that there might be overlap between components.}

\item{spls_selection_method}{Indicates the method for variable selection on sPLS results. One of 'correlation' or 'loadings' (default: 'loadings'). Loadings method will choose those variables maximizing variance across the samples, while correlation method will choose those variables with a higher correlation with other variables, that is those variables more distant to the origin in the correlation plot.}

\item{GO_similarity_measure}{The similarity measured employed to cluster and visualize GO terms, one in 'Resnik', 'Lin', 'Rel', 'Jiang', 'Wang', 'UI', 'binary', 'bray-curtis', 'cosine', 'all' (default: 'Wang')}

\item{GOA_search_universe}{The GOA universe to use in functional similarity measures (i.e.: only aplicable to 'binary', 'bray-curtis', 'cosine', 'UI'), one in 'human', 'uniprot', 'gene_list' (default: 'human')}

\item{enrichment_significance_thr}{The threshold to consider a GO enrichment result as significant (default: 0.01)}

\item{GO_ontology}{The GO ontology on which to perform the enrichment, one of 'BP', 'MF', 'CC' (default: BP)}

\item{multiple_test_adjustment}{Flag to run multiple test adjustment on GO enrichment results (default: F)}
}
\description{
Runs the TCGAome analysis pipeline on the given tumor types.
}
\examples{
run_TCGAome(c('BRCA', 'OV'))
}
\keyword{TCGAome}

