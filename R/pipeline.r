
results_evaluation <- function(spls_results, mcia_results)
{
  flog.info("Evaluating results from MCIA and sPLS ...")

  output_dir = paste(RESULTS_FOLDER, "results_evaluation", sep="/")
  dir.create(output_dir)

  # Compares to TSGs and oncogenes from publication Vogelstein et al. (2013). Cancer genome landscapes. Science
  vogelstein_oncogenes_and_tsg = read.csv(paste(get_package_folder("inst"), "vogelstein_oncogenes_and_tsg.processed.txt", sep="/"), sep="\t")
  #oncogenes = as.character(oncogenes_and_tsg[oncogenes_and_tsg$Classification == "Oncogene", 1])
  #tsgs = as.character(oncogenes_and_tsg[oncogenes_and_tsg$Classification == "TSG", 1])

  # Plots Venn diagram
  plot_triple_venn(
    mcia_results$selected_variables,
    spls_results$selected_variables,
    as.character(vogelstein_oncogenes_and_tsg$Gene.Symbol),
    c("MCIA selected genes", "sPLS selected genes", "Oncogenes and TSGs (Vogelstein et al. 2013)"),
    output_dir=output_dir,
    file_name="selected_genes_vs_vogelstein_venn.png"
  )

  # Compares to Tumor Suppressor Gene Database selected by differential gene expression pan-cancer analysis
  tsg_database = read.table(paste(get_package_folder("inst"), "tsg_database_sig_exp.txt", sep="/"), header=T, sep="\t")

  # Plots Venn diagram
  plot_triple_venn(
    mcia_results$selected_variables,
    spls_results$selected_variables,
    tsg_database$Symbol,
    c("MCIA selected genes", "sPLS selected genes", "Tumor Suppressor Genes Database"),
    output_dir=output_dir,
    file_name="selected_genes_vs_tsgdb_venn.png"
  )


  intersect_mcia_vogelstein = intersect(vogelstein_oncogenes_and_tsg$Gene.Symbol, mcia_results$selected_variables)
  intersect_spls_vogelstein = intersect(vogelstein_oncogenes_and_tsg$Gene.Symbol, spls_results$selected_variables)
  intersect_mcia_tsgdb = intersect(tsg_database$Symbol, mcia_results$selected_variables)
  intersect_spls_tsgdb = intersect(tsg_database$Symbol, spls_results$selected_variables)
  gene_universe = gsub("\\.df.", "", row.names(mcia_results$results$mcoa$axis))
  intersect_vogelstein_geneuniverse = intersect(vogelstein_oncogenes_and_tsg$Gene.Symbol, gene_universe)
  intersect_tsgdb_geneuniverse = intersect(tsg_database$Symbol, gene_universe)

  # Calculates the probability of finding a TSG or oncogene in the selected genes by chance
  trials_mcia = length(mcia_results$selected_variables)
  trials_spls = length(spls_results$selected_variables)
  n = length(gene_universe)
  n1 = length(intersect_vogelstein_geneuniverse)
  n2 = length(intersect_tsgdb_geneuniverse)

  binomial_coefficient <- function(n, k) exp(lfactorial(n) - lfactorial(k) - lfactorial(n-k))

  # Probability of finding 1 of the selected genes being part of the Vogelstein list
  p1_mcia = 1 - exp(log(binomial_coefficient(n-n1, trials_mcia)) - log(binomial_coefficient(n, trials_mcia)))
  flog.info("Probability that 1 of the MCIA selected variables is part of the Vogelstein list by chance: %f", p1_mcia)
  p1_spls = 1 - exp(log(binomial_coefficient(n-n1, trials_spls)) - log(binomial_coefficient(n, trials_spls)))
  flog.info("Probability that 1 of the sPLS selected variables is part of the Vogelstein list by chance: %f", p1_spls)

  # Probability of finding 1 of the selected genes being part of the TSG database
  p2_mcia = 1 - exp(log(binomial_coefficient(n-n2, trials_mcia)) - log(binomial_coefficient(n, trials_mcia)))
  flog.info("Probability that 1 of the MCIA selected variables is part of the Tumor Suppresor Gene Database by chance: %f", p2_mcia)
  p2_spls = 1 - exp(log(binomial_coefficient(n-n2, trials_spls)) - log(binomial_coefficient(n, trials_spls)))
  flog.info("Probability that 1 of the sPLS selected variables is part of the Tumor Suppresor Gene Database by chance: %f", p2_spls)

  # Approximation of the probability of finding N of the selected genes being part of the Vogelstein list
  p1_n_mcia = p1_mcia^length(intersect_mcia_vogelstein)
  flog.info("MCIA selected genes contains %d genes in the Vogelstein list", length(intersect_mcia_vogelstein))
  flog.info("Probability that %d of the MCIA selected variables is part of the Vogelstein list by chance: %f", length(intersect_mcia_vogelstein), p1_n_mcia)
  p1_n_spls = p1_spls^length(intersect_spls_vogelstein)
  flog.info("sPLS selected genes contains %d genes in the Vogelstein list", length(intersect_spls_vogelstein))
  flog.info("Probability that %d of the sPLS selected variables is part of the Vogelstein list by chance: %f", length(intersect_spls_vogelstein), p1_n_spls)

  # Approximation of the probability of finding N of the selected genes being part of the TSG database
  p2_n_mcia = p2_mcia^length(intersect_mcia_tsgdb)
  flog.info("MCIA selected genes contains %d genes in the Tumor Suppressor Gene Database", length(intersect_mcia_tsgdb))
  flog.info("Probability that %d of the MCIA selected variables is part of the Tumor Suppressor Gene Database by chance: %f", length(intersect_mcia_tsgdb), p2_n_mcia)
  p2_n_spls = p2_spls^length(intersect_spls_tsgdb)
  flog.info("sPLS selected genes contains %d genes in the Tumor Suppressor Gene Database", length(intersect_spls_tsgdb))
  flog.info("Probability that %d of the sPLS selected variables is part of the Tumor Suppressor Gene Database by chance: %f", length(intersect_spls_tsgdb), p2_n_spls)

  attributes = c("Total number of genes evaluated",
                 "No. of selected variables by MCIA",
                 "No. of selected variables by sPLS",
                 "No. of genes in the Vogelstein list",
                 "No. of genes in the Tumor Suppressor Gene Database",
                 "No. of MCIA matches in the Vogelstein list",
                 "No. of sPLS matches in the Vogelstein list",
                 "No. of MCIA matches in the Tumor Suppressor Gene Database",
                 "No. of sPLS matches in the Tumor Suppressor Gene Database",
                 "Probability of finding 1 match on MCIA selected genes on the Vogelstein list",
                 "Probability of finding 1 match on sPLS selected genes on the Vogelstein list",
                 "Probability of finding 1 match on MCIA selected genes on the Tumor Suppressor Gene Database",
                 "Probability of finding 1 match on sPLS selected genes on the Tumor Suppressor Gene Database",
                 "Probability of finding N matches on MCIA selected genes on the Vogelstein list",
                 "Probability of finding N matches on sPLS selected genes on the Vogelstein list",
                 "Probability of finding N matches on MCIA selected genes on the Tumor Suppressor Gene Database",
                 "Probability of finding N matches on sPLS selected genes on the Tumor Suppressor Gene Database"
                 )
  values = c(length(gene_universe),
             length(mcia_results$selected_variables),
             length(spls_results$selected_variables),
             length(intersect_vogelstein_geneuniverse),
             length(intersect_tsgdb_geneuniverse),
             length(intersect_mcia_vogelstein),
             length(intersect_spls_vogelstein),
             length(intersect_mcia_tsgdb),
             length(intersect_spls_tsgdb),
             p1_mcia,
             p1_spls,
             p2_mcia,
             p2_spls,
             p1_n_mcia,
             p1_n_spls,
             p2_n_mcia,
             p2_n_spls
             )

  write.table(data.frame(attribute=attributes, value=values), file = paste(output_dir, "results_evaluation.txt", sep="/"), quote=F, row.names=F, sep = "\t")

  flog.info("Results evaluation finished.")
}


#' Returns the supported tumor types.
#'
#' @keywords TCGAome
#' @export
#' @examples
#' get_TCGAome_tumor_types()
get_TCGAome_tumor_types <- function(){
  names(RPPA_ANNOTATIONS)
}


#' Runs the TCGAome analysis pipeline on the given tumor types.
#' @param tumor_types Vector of tumor types to be analyzed from those available at TCGA (run RTCGAToolbox::getFirehoseDatasets() to see all available types)
#' @param run_pca_analysis Flag indicating if the PCA should run, this is part of the data pre-analysis (default: TRUE)
#' @param run_hclust_analysis Flag indicating if the hierarchical clustering analysis should run, this is part of the data pre-analysis (default: FALSE)
#' @param run_rgcca Flag indicating if the Regularized Generalized Canonical Correlation Analysis (RGCCA) should run (default: FALSE). Beware that this analysis is not intended for datasets with number of variables >> number of samples.
#' @param run_rcca Flag indicating if the Regularized Canonical Correlation Analysis (rCCA) should run (default: FALSE). Beware that this analysis is not intended for datasets with number of variables >> number of samples.
#' @param topN Indicates the top number of variables to select on MCIA and sPLS results (default: 5). It will select N variables on each of the data types, on each of the three first components and on each extreme of range, that is a maximum of 2*3*2*N, considering that there might be overlap between components.
#' @param spls_selection_method Indicates the method for variable selection on sPLS results. One of "correlation" or "loadings" (default: "loadings"). Loadings method will choose those variables maximizing variance across the samples, while correlation method will choose those variables with a higher correlation with other variables, that is those variables more distant to the origin in the correlation plot.
#' @param GO_similarity_measure The similarity measured employed to cluster and visualize GO terms, one in "Resnik", "Lin", "Rel", "Jiang", "Wang", "UI", "binary", "bray-curtis", "cosine", "all" (default: "Wang")
#' @param GOA_search_universe The GOA universe to use in functional similarity measures (i.e.: only aplicable to "binary", "bray-curtis", "cosine", "UI"), one in "human", "uniprot", "gene_list" (default: "human")
#' @param enrichment_significance_threshold The threshold to consider a GO enrichment result as significant (default: 0.01)
#' @param GO_ontology The GO ontology on which to perform the enrichment, one of "BP", "MF", "CC" (default: BP)
#' @param multiple_test_adjustment Flag to run multiple test adjustment on GO enrichment results (default: F)
#'
#' @keywords TCGAome
#' @export
#' @examples
#' run_TCGAome(c("BRCA", "OV"))
run_TCGAome <- function(tumor_types,
                        run_pca_analysis=TRUE,
                        run_hclust_analysis=FALSE,
                        run_rgcca=FALSE,
                        run_rcca=FALSE,
                        topN = 5,
                        spls_selection_method = "loadings",
                        GO_similarity_measure = "Wang",
                        GOA_search_universe = "human",
                        enrichment_significance_threshold = 0.01,
                        GO_ontology = "BP",
                        multiple_test_adjustment = FALSE
                        ){

  # Creates and sets the results folder for the current run
  RESULTS_FOLDER <<- get_results_folder()

  # Configures logging
  configure_logging(RESULTS_FOLDER)

  flog.info("Results directory: %s", RESULTS_FOLDER)

  # Downloads data for RNAseq and RPPA
  matrices = get_tcga_data(tumor_types)

  # Preprocessing
  preprocessed_matrices = preprocess_data(X=matrices$X, Y=matrices$Y, correlation.thr = 0.7)

  # Pre-analysis
  descriptive_analysis(X=preprocessed_matrices$X, Y=preprocessed_matrices$Y, Z=matrices$Z)

  if (run_pca_analysis){
    pca_analysis(X=preprocessed_matrices$X, Y=preprocessed_matrices$Y, Z=matrices$Z)
  } else {
    flog.info("Principal Component Analysis disabled.")
  }

  if (run_hclust_analysis){
    hclust_analysis(X=preprocessed_matrices$X, Y=preprocessed_matrices$Y, Z=matrices$Z)
  } else {
    flog.info("Hierarchichal clustering analysis disabled.")
  }

  # Runs MCIA
  mcia_results = mcia_analysis(X = preprocessed_matrices$X, Y = preprocessed_matrices$Y, Z = matrices$Z, topN = topN, cia.nf = 5)

  # Runs sPLS
  spls_results = spls_analysis(X = preprocessed_matrices$X, Y = preprocessed_matrices$Y, Z = matrices$Z, topN = topN, selection_method = spls_selection_method)


  # Runs RGCCA
  if (run_rgcca){
    library(RGCCA)
    rgcca_results = rgcca_analysis(X=preprocessed_matrices$X, Y=preprocessed_matrices$Y, Z=matrices$Z)
  } else {
    flog.info("Regularized Generalized Canonical Correlation Analysis disabled.")
  }

  # Runs rCCA
  #preprocessed_matrices = preprocess.data(matrices$X, matrices$Y, correlation.thr = 0.5)
  if (run_rcca){
    rcaa_results = rcaa_analysis(preprocessed_matrices$X, preprocessed_matrices$Y, Z=matrices$Z)
  } else {
    flog.info("Regularized Canonical Correlation Analysis disabled.")
  }

  # Evaluates results
  results_evaluation(mcia_results = mcia_results, spls_results = spls_results)

  # Performs pathway enrichment analysis
  functional_analysis(mcia_results = mcia_results,
                      spls_results = spls_results,
                      GO_similarity_measure = GO_similarity_measure,
                      GOA_search_universe = GOA_search_universe,
                      enrichment_significance_threshold = enrichment_significance_threshold,
                      GO_ontology = GO_ontology,
                      multiple_test_adjustment = multiple_test_adjustment)


  # TODO: enrich list of known genes and compare to results from MCIA and sPLS

  ###
  ## There is no way to run this analysis we need to reduce dimensions
  ###
  # TODO: create a test mode that reduces matrices
  # Random selection of variables (10 genes and 10 proteins)
  # X.random.subset = X[, round(runif(10, min = 1, max=dim(X)[2]))]
  # Y.random.subset = Y[, round(runif(10, min = 1, max=dim(Y)[2]))]

  # Reads a list of Tumor Supressor Genes for BRCA tumor type
  #brca_tsg = read.table("./data/BRCA_down_regulated_TSgenes.txt", col.names = c("gene_id", "gene_name"))
  #brca_tsg = intersect(colnames(X), brca_tsg$gene_name)
  #X.tsg.subset = X[, brca_tsg]
  #X.random.subset = X[, brca_tsg[round(runif(100, min = 1, max=length(brca_tsg)))]]
  #Y.random.subset = Y[, round(runif(10, min = 1, max=dim(Y)[2]))]


  # Create .md, .html, and .pdf files
  #knit("File.Rmd")
  #markdownToHTML('File.md', 'File.html', options=c("use_xhml"))
  #system("pandoc -s File.html -o File.pdf")


  flog.info("TCGAome finished!")
  flog.info("Your results are at %s", RESULTS_FOLDER)
}



