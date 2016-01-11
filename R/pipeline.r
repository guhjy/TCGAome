#library(ReactomePA)
#library(RTCGAToolbox)
#library(mixOmics)
#library(omicade4)
#library(RGCCA)
#library(dplyr)
#library(VennDiagram)







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


###
### Functional analysis
###
functional_analysis <- function (spls_results, mcia_results)
{

  # Loads ensembl biomart for Homo sapiens


  oncogenes_and_tsg = read.csv("./data/oncogenes_and_tsg.processed.txt", sep="\t")

  # get entrezid gene list from genes (it loses one gene on conversion, I guess it is "C11orf75")
  mcia_selected_variables_entrezid = query_biomart(attributes=c("entrezgene"), filters=c("hgnc_symbol"), values=mcia_results$selected_variables)[, 1]
  spls_selected_variables_entrezid = query_biomart(attributes=c("entrezgene"), filters=c("hgnc_symbol"), values=spls_results$selected_variables)[, 1]
  oncogenes_and_tsg_entrezid = query_biomart(attributes=c("entrezgene"), filters=c("hgnc_symbol"), values=oncogenes_and_tsg$Gene.Symbol)[, 1]

  # Compute enrichment analysis on Reactome
  mcia_enrichment = enrichPathway(gene=mcia_selected_variables_entrezid, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T)
  spls_enrichment = enrichPathway(gene=spls_selected_variables_entrezid, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T)
  oncogenes_and_tsg_enrichment = enrichPathway(gene=oncogenes_and_tsg_entrezid, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T)

  # Plots Venn diagram
  plot_triple_venn(
    mcia_enrichment@result$ID,
    spls_enrichment@result$ID,
    oncogenes_and_tsg_enrichment@result$ID,
    c("MCIA", "sPLS", "Oncogenes and TSGs"),
    "results", "triple_venn.png"
  )

  # Calculates GO enrichment for MCIA
  ENRICHMENT_SIGNIFICANCE_THRESHOLD = 0.05
  ONTOLOGY = "BP"
  mcia_go_enrichment_results = get_go_enrichment(mcia_selected_variables_entrezid,
                                                 pvalue_threshold = ENRICHMENT_SIGNIFICANCE_THRESHOLD,
                                                 ontology = ONTOLOGY)
  mcia_go_enrichment = mcia_go_enrichment_results$data.frame
  mcia_TopGOdata = mcia_go_enrichment_results$TopGOdata


  # Calculates frequency (size) of GO terms according to GOA (we use the whole uniprot to calculate size)
  mcia_go_enrichment$size = as.vector(sapply(mcia_go_enrichment$GO, FUN = get_goa_size))


  write.table(mcia_go_enrichment, file = "results/GO_enrichment/MCIA/go_enrichment.txt", sep="\t", row.names = F, quote = F)

  #TODO: separate by GO subontology (i.e.: biological process, etc.)
  #mcia_go_details = query_biomart(attributes=c("go_id", "name_1006", "definition_1006", "namespace_1003"), filters=c("go_id"), values=mcia_go_enrichment$GO)
  # join
  # filter out other than bp

  # Plots enrichment for MCIA results

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/Wang",
                   method = "Wang",
                   clustering_method = "pam",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/binary",
                   method = "binary",
                   clustering_method = "pam",
                   search_universe = "human",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/binary_gene_list",
                   method = "binary",
                   clustering_method = "pam",
                   search_universe = "gene_list",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/binary_uniprot",
                   method = "binary",
                   clustering_method = "pam",
                   search_universe = "uniprot",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/UI",
                   method = "UI",
                   clustering_method = "pam",
                   search_universe = "human",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/UI_gene_list",
                   method = "UI",
                   clustering_method = "pam",
                   search_universe = "gene_list",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/UI_uniprot",
                   method = "UI",
                   clustering_method = "pam",
                   search_universe = "uniprot",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/bray-curtis",
                   method = "bray-curtis",
                   clustering_method = "pam",
                   search_universe = "human",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/bray-curtis_gene_list",
                   method = "bray-curtis",
                   clustering_method = "pam",
                   search_universe = "gene_list",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/bray-curtis_uniprot",
                   method = "bray-curtis",
                   clustering_method = "pam",
                   search_universe = "uniprot",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/cosine",
                   method = "cosine",
                   clustering_method = "pam",
                   search_universe = "human",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/cosine_gene_list",
                   method = "cosine",
                   clustering_method = "pam",
                   search_universe = "gene_list",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/cosine_uniprot",
                   method = "cosine",
                   clustering_method = "pam",
                   search_universe = "uniprot",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/Resnik",
                   method = "Resnik",
                   clustering_method = "pam",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/Lin",
                   method = "Lin",
                   clustering_method = "pam",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/Rel",
                   method = "Rel",
                   clustering_method = "pam",
                   ont=ONTOLOGY)

  cluster_and_plot(enrichment_results = mcia_go_enrichment,
                   gene_list = mcia_results$selected_variables,
                   TopGOdata = mcia_TopGOdata,
                   output_dir="results/GO_enrichment/MCIA/Jiang",
                   method = "Jiang",
                   clustering_method = "pam",
                   ont=ONTOLOGY)



  # Retrieve Reactome and biological process GO terms for MCIA results
  #mcia_reactome = query_biomart(attributes=c("entrezgene", "hgnc_symbol", "reactome"), filters=c("entrezgene"), values=mcia_selected_variables_entrezid)
  #mcia_go = query_biomart(attributes=c("entrezgene", "hgnc_symbol", "go_id", "name_1006", "definition_1006", "namespace_1003"), filters=c("entrezgene"), values=mcia_selected_variables_entrezid)
  #mcia_go_id_unique = unique(mcia_go$go_id)

  # Keeps only the most specific terms, those not having descendents in the list
  most_specific_terms = get_most_specific_terms(mcia_go_enrichment$GO)
  mcia_go_enrichment = mcia_go_enrichment[mcia_go_enrichment$GO %in% most_specific_terms, ]

  # Counts ocurrences of each term
  mcia_go_bp = mcia_go[mcia_go$namespace_1003 == 'biological_process' & mcia_go$go_id %in% most_specific_terms, ]
  mcia_go_bp_count = as.data.frame(table(mcia_go_bp$go_id))
  mcia_go_bp_count = mcia_go_bp_count[order(-mcia_go_bp_count$Freq), ]

  # Top 50 ????
  mcia_go_bp_top50 = mcia_go_bp_count[1:50, ]

  # plots frequency distribution of top 50 GO biological processes
  barplot(mcia_go_bp_top50$Freq, xlab="GO biological process", ylab="frequency", names.arg = mcia_go_bp_top50$Var1)
  mcia_go_bp_top50$Var1 = factor(mcia_go_bp_top50$Var1,
                                 levels = mcia_go_bp_top50$Var1[order(mcia_go_bp_top50$Freq, decreasing=T)])
  q <- qplot(x=mcia_go_bp_top50$Var1, y=mcia_go_bp_top50$Freq,
             data=mcia_go_bp_top50, geom="bar", stat="identity", xlab = "GO biological process", ylab="No. of genes")
  q + theme(axis.text.x = element_text(angle = 90, vjust=0, hjust = 1))


  # TODO: Retrieve Reactome and GO terms for sPLS results
  spls_pathways = getBM(attributes=c("entrezgene", "hgnc_symbol", "reactome", "go_id", "name_1006", "definition_1006"), filters=c("entrezgene"), values=spls_selected_variables_entrezid, mart=ensembl)


}


# Evaluates the similarity measure and clustering results
go_clustering_evaluation <- function (method){

  base_folder = "results/GO_enrichment"

  # Retrieves results for distance to centroid of each similarity metric
  results_binary = read_clustering_results(paste(method, "binary", sep="/"), "binary")
  dist2centroid = results_binary$dist2centroid

  results_UI = read_clustering_results(paste(method, "UI", sep="/"), "UI")
  dist2centroid = merge(dist2centroid, results_UI$dist2centroid)

  results_BC = read_clustering_results(paste(method, "bray-curtis", sep="/"), "BC")
  dist2centroid = merge(dist2centroid, results_BC$dist2centroid)

  results_cosine = read_clustering_results(paste(method, "cosine", sep="/"), "cosine")
  dist2centroid = merge(dist2centroid, results_cosine$dist2centroid)

  results_binary_intra = read_clustering_results(paste(method, "binary_gene_list", sep="/"), "binary_intra")
  dist2centroid = merge(dist2centroid, results_binary_intra$dist2centroid)

  results_UI_intra = read_clustering_results(paste(method, "UI_gene_list", sep="/"), "UI_intra")
  dist2centroid = merge(dist2centroid, results_UI_intra$dist2centroid)

  results_BC_intra = read_clustering_results(paste(method, "bray-curtis_gene_list", sep="/"), "BC_intra")
  dist2centroid = merge(dist2centroid, results_BC_intra$dist2centroid)

  results_cosine_intra = read_clustering_results(paste(method, "cosine_gene_list", sep="/"), "cosine_intra")
  dist2centroid = merge(dist2centroid, results_cosine_intra$dist2centroid)

  results_jiang = read_clustering_results(paste(method, "Jiang", sep="/"), "Jiang")
  dist2centroid = merge(dist2centroid, results_jiang$dist2centroid)

  results_lin = read_clustering_results(paste(method, "Lin", sep="/"), "Lin")
  dist2centroid = merge(dist2centroid, results_lin$dist2centroid)

  results_rel = read_clustering_results(paste(method, "Rel", sep="/"), "Schliker")
  dist2centroid = merge(dist2centroid, results_rel$dist2centroid)

  results_resnik = read_clustering_results(paste(method, "Resnik", sep="/"), "Resnik")
  dist2centroid = merge(dist2centroid, results_resnik$dist2centroid)

  results_wang = read_clustering_results(paste(method, "Wang", sep="/"), "Wang")
  dist2centroid = merge(dist2centroid, results_wang$dist2centroid)


  # Plots boxplot
  dist2centroid_pivot = melt(dist2centroid, id.vars='GO', measure.vars=names(dist2centroid)[2:14])
  png(paste(paste(base_folder, method, sep="/"),"boxplot_dist2centroid.png", sep="/"))
  ggplot(dist2centroid_pivot, aes(x=variable, y=value)) +
    geom_boxplot() +
    xlab("Metric") +
    ylab("Distance to centroid") +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5)) +
    stat_summary(fun.y = mean, geom="point",colour="darkred", size=2)
  dev.off()


  # Retrieves clustering results metrics for all similarity measures
  clustering_metrics = as.data.frame(bind_rows(list(
    get_clustering_metrics("binary", results_binary$data),
    get_clustering_metrics("UI", results_UI$data),
    get_clustering_metrics("BC", results_BC$data),
    get_clustering_metrics("cosine", results_cosine$data),
    get_clustering_metrics("binary_intra", results_binary_intra$data),
    get_clustering_metrics("UI_intra", results_UI_intra$data),
    get_clustering_metrics("BC_intra", results_BC_intra$data),
    get_clustering_metrics("cosine_intra", results_cosine_intra$data),
    get_clustering_metrics("Jiang", results_jiang$data),
    get_clustering_metrics("Lin", results_lin$data),
    get_clustering_metrics("Rel", results_rel$data),
    get_clustering_metrics("Resnik", results_resnik$data),
    get_clustering_metrics("Wang", results_wang$data)
  )))


  # Plots barplot with number of clusters
  png(paste(paste(base_folder, method, sep="/"),"count_clusters.png", sep="/"))
  ggplot(clustering_metrics, aes(x = factor(measure, levels=measure), y = clusters)) +
    geom_bar(stat = "identity", position=position_dodge(width = 0.9), width=0.5) +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5)) +
    ylab("# of clusters") +
    xlab("Measures")
  dev.off()

  clustering_metrics_pivot = melt(clustering_metrics, id.vars='measure', measure.vars=names(clustering_metrics)[4:8])
  clustering_metrics_pivot$value = as.numeric(clustering_metrics_pivot$value)
  clustering_metrics_pivot$measure = factor(clustering_metrics_pivot$measure, levels=clustering_metrics_pivot$measure)

  png(paste(paste(base_folder, method, sep="/"),"descriptive_stats_clusters.png", sep="/"))
  ggplot(clustering_metrics_pivot, aes(x = measure, y = value, fill=factor(variable))) +
    geom_bar(stat = "identity", position=position_dodge(width = 0.9), width=0.5) +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5)) +
    ylab("Value") +
    xlab("Measures") +
    labs(fill="")
  #coord_flip() +
  #scale_x_discrete(limits = rev(levels(clustering_metrics_pivot$measure)))
  #scale_fill_continuous(low = "grey", high = "red", space = "Lab", name = "g = 0")
  dev.off()



  wang_set = unique(results_wang$data$cluster)
  resnik_set = unique(results_resnik$data$cluster)

  cosine_set = unique(results_cosine$data$cluster)
  png(paste(paste(base_folder, method, sep="/"),"venn_cosine_wang_resnik.png", sep="/"),width=800,height=700, res=96)
  plot_venn(wang_set, resnik_set, cosine_set, c("Wang", "Resnik", "Cosine"))
  dev.off()

  UI_set = unique(results_UI$data$cluster)
  png(paste(paste(base_folder, method, sep="/"),"venn_UI_wang_resnik.png", sep="/"),width=800,height=700, res=96)
  plot_venn(wang_set, resnik_set, UI_set, c("Wang", "Resnik", "UI"))
  dev.off()

  binary_set = unique(results_binary$data$cluster)
  png(paste(paste(base_folder, method, sep="/"),"venn_binary_wang_resnik.png", sep="/"),width=800,height=700, res=96)
  plot_venn(wang_set, resnik_set, binary_set, c("Wang", "Resnik", "binary"))
  dev.off()

  BC_set = unique(results_BC$data$cluster)
  png(paste(paste(base_folder, method, sep="/"),"venn_BC_wang_resnik.png", sep="/"),width=800,height=700, res=96)
  plot_venn(wang_set, resnik_set, BC_set, c("Wang", "Resnik", "BC"))
  dev.off()

  cosine_intra_set = unique(results_cosine_intra$data$cluster)
  png(paste(paste(base_folder, method, sep="/"),"venn_cosineintra_wang_resnik.png", sep="/"),width=800,height=700, res=96)
  plot_venn(wang_set, resnik_set, cosine_intra_set, c("Wang", "Resnik", "Cosine intra"))
  dev.off()

  UI_intra_set = unique(results_UI_intra$data$cluster)
  png(paste(paste(base_folder, method, sep="/"),"venn_UIintra_wang_resnik.png", sep="/"),width=800,height=700, res=96)
  plot_venn(wang_set, resnik_set, UI_intra_set, c("Wang", "Resnik", "UI intra"))
  dev.off()

  binary_intra_set = unique(results_binary_intra$data$cluster)
  png(paste(paste(base_folder, method, sep="/"),"venn_binaryintra_wang_resnik.png", sep="/"),width=800,height=700, res=96)
  plot_venn(wang_set, resnik_set, binary_intra_set, c("Wang", "Resnik", "binary intra"))
  dev.off()

  BC_intra_set = unique(results_BC_intra$data$cluster)
  png(paste(paste(base_folder, method, sep="/"),"venn_BCintra_wang_resnik.png", sep="/"),width=800,height=700, res=96)
  plot_venn(wang_set, resnik_set, BC_intra_set, c("Wang", "Resnik", "BC intra"))
  dev.off()
}


# Reads cluster results
read_clustering_results <- function(folder, measure){

  base_folder = "results/GO_enrichment"
  folder = paste(base_folder, folder, sep="/")
  dist2centroid_file = paste(folder, "distances_to_centroid.txt", sep="/")
  dist2centroid = read.table(dist2centroid_file, header = T, sep = "\t", stringsAsFactors = F)
  names(dist2centroid) = c("GO", measure)
  results_file = paste(folder, "enrichment_results.txt", sep="/")
  data = read.table(results_file, header = T, sep = "\t", stringsAsFactors = F)
  list(dist2centroid=dist2centroid, data=data)
}

# Calculates descriptive stats on clustering results
get_clustering_metrics <- function(measure, data){
  clustering_summary = table(data$cluster)
  data.frame(
    measure=measure,
    nodes=length(data$GO),
    clusters=length(unique(data$cluster)),
    mean=round(mean(clustering_summary)),
    median=round(median(clustering_summary)),
    mode=names(which.max(table(clustering_summary))),
    min=min(clustering_summary),
    max=max(clustering_summary)
  )
}


plot_venn <- function(set1, set2, set3, labels){


  print({
    draw.triple.venn(area1 = length(set1),
                     area2 = length(set2),
                     area3 = length(set3),
                     n12 = length(intersect(set1, set2)),
                     n23 = length(intersect(set2, set3)),
                     n13 = length(intersect(set1, set3)),
                     n123 = length(intersect(set1, intersect(set2, set3))),
                     category = labels, lty = "blank",
                     fill = c(red500, green500, indigo500))
  })
}


#' Runs the TCGAome analysis pipeline on the given tumor types.
#' @param tumor_types Vector of tumor types to be analyzed from those available at TCGA (run RTCGAToolbox::getFirehoseDatasets() to see all available types)
#' @param run_pca_analysis Flag indicating if the PCA should run, this is part of the data pre-analysis (default: TRUE)
#' @param run_hclust_analysis Flag indicating if the hierarchical clustering analysis should run, this is part of the data pre-analysis (default: FALSE)
#' @param run_rgcca Flag indicating if the Regularized Generalized Canonical Correlation Analysis (RGCCA) should run (default: FALSE). Beware that this analysis is not intended for datasets with number of variables >> number of samples.
#' @param run_rcca Flag indicating if the Regularized Canonical Correlation Analysis (rCCA) should run (default: FALSE). Beware that this analysis is not intended for datasets with number of variables >> number of samples.
#' @param topN Indicates the top number of variables to select on MCIA and sPLS results (default: 5). It will select N variables on each of the data types, on each of the three first components and on each extreme of range, that is a maximum of 2*3*2*N, considering that there might be overlap between components.
#' @param spls_selection_method Indicates the method for variable selection on sPLS results. One of "correlation" or "loadings" (default: "loadings"). Loadings method will choose those variables maximizing variance across the samples, while correlation method will choose those variables with a higher correlation with other variables, that is those variables more distant to the origin in the correlation plot.
#'
#' @keywords TCGAome
#' @export
#' @examples
#' run.TCGAome(c("BRCA", "OV"))
run.TCGAome <- function(tumor_types,
                        run_pca_analysis=TRUE,
                        run_hclust_analysis=FALSE,
                        run_rgcca=FALSE,
                        run_rcca=FALSE,
                        variable_selection_topN = 5,
                        spls_selection_method = "loadings"
                        ){

  # Loads bioconductor packages
  loads_dependencies()

  # Creates and sets the results folder for the current run
  RESULTS_FOLDER <<- get_results_folder()

  # Configures logging
  configure_logging(RESULTS_FOLDER)

  flog.info("Results directory: %s", RESULTS_FOLDER)

  # Downloads data for RNAseq and RPPA
  #matrices = get.data(c("BRCA", "OV"))
  # tumor_types = c("BRCA", "OV")
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
  functional_analysis(mcia_results = mcia_results, spls_results = spls_results)


  ###
  ## There is no way to run this analysis we need to reduce dimensions
  ###

  # Random selection of variables (10 genes and 10 proteins)
  # X.random.subset = X[, round(runif(10, min = 1, max=dim(X)[2]))]
  # Y.random.subset = Y[, round(runif(10, min = 1, max=dim(Y)[2]))]

  # Reads a list of Tumor Supressor Genes for BRCA tumor type
  brca_tsg = read.table("./data/BRCA_down_regulated_TSgenes.txt", col.names = c("gene_id", "gene_name"))
  brca_tsg = intersect(colnames(X), brca_tsg$gene_name)
  X.tsg.subset = X[, brca_tsg]
  X.random.subset = X[, brca_tsg[round(runif(100, min = 1, max=length(brca_tsg)))]]
  Y.random.subset = Y[, round(runif(10, min = 1, max=dim(Y)[2]))]


  # Create .md, .html, and .pdf files
  knit("File.Rmd")
  markdownToHTML('File.md', 'File.html', options=c("use_xhml"))
  system("pandoc -s File.html -o File.pdf")
}



