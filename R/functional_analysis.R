###
### Functional analysis
###
functional_analysis <- function (spls_results, mcia_results,
                                 GO_similarity_measure = "Wang",
                                 GOA_search_universe = "human",
                                 enrichment_significance_threshold = 0.05,
                                 GO_ontology = "BP",
                                 multiple_test_adjustment = FALSE)
{

  #oncogenes_and_tsg = read.csv("./data/oncogenes_and_tsg.processed.txt", sep="\t")

  # get entrezid gene list from genes (it loses one gene on conversion, I guess it is "C11orf75")
  mcia_selected_variables_entrezid = query_biomart(attributes=c("entrezgene"), filters=c("hgnc_symbol"), values=mcia_results$selected_variables)[, 1]
  spls_selected_variables_entrezid = query_biomart(attributes=c("entrezgene"), filters=c("hgnc_symbol"), values=spls_results$selected_variables)[, 1]
  #oncogenes_and_tsg_entrezid = query_biomart(attributes=c("entrezgene"), filters=c("hgnc_symbol"), values=oncogenes_and_tsg$Gene.Symbol)[, 1]

  # Compute enrichment analysis on Reactome
  #mcia_enrichment = enrichPathway(gene=mcia_selected_variables_entrezid, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T)
  #spls_enrichment = enrichPathway(gene=spls_selected_variables_entrezid, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T)
  #oncogenes_and_tsg_enrichment = enrichPathway(gene=oncogenes_and_tsg_entrezid, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T)

  # Plots Venn diagram
  #plot_triple_venn(
  #  mcia_enrichment@result$ID,
  #  spls_enrichment@result$ID,
  #  oncogenes_and_tsg_enrichment@result$ID,
  #  c("MCIA", "sPLS", "Oncogenes and TSGs"),
  #  "results", "triple_venn.png"
  #)

  flog.info("Functional analysis starting ...")

  flog.info("... on MCIA results ...")

  # Runs functional analysis on MCIA results
  output_dir = paste(RESULTS_FOLDER, "MCIA/functional_analysis", sep="/")
  dir.create(output_dir, recursive=T)

  # Calculates GO enrichment for MCIA
  mcia_go_enrichment_results = get_go_enrichment(gene_list = mcia_selected_variables_entrezid,
                                                 pvalue_threshold = enrichment_significance_threshold,
                                                 ontology = GO_ontology,
                                                 multiple_test_adjustment = multiple_test_adjustment)
  mcia_go_enrichment = mcia_go_enrichment_results$data.frame
  mcia_TopGOdata = mcia_go_enrichment_results$TopGOdata

  if (dim(mcia_go_enrichment)[1] > 0) {
    write.table(mcia_go_enrichment, file = paste(output_dir, "go_enrichment.txt", sep="/"), sep="\t", row.names = F, quote = F)
    flog.info("GO enrichment on MCIA selected genes returned %d enriched GO terms.", dim(mcia_go_enrichment)[1])
    # Plots enrichment for MCIA results
    cluster_and_plot(enrichment_results = mcia_go_enrichment,
                     gene_list = mcia_results$selected_variables,
                     TopGOdata = mcia_TopGOdata,
                     output_dir = output_dir,
                     method = GO_similarity_measure,
                     clustering_method = "pam",
                     ont = GO_ontology,
                     search_universe = GOA_search_universe)
  } else {
    flog.info("GO enrichment on MCIA selected genes gave no results.")
    if (multiple_test_adjustment) {
      flog.info("Try disabling multiple test adjustment.")
    }
  }



  flog.info("... on sPLS results ...")

  # Runs functional analysis on sPLS results
  output_dir = paste(RESULTS_FOLDER, "sPLS/functional_analysis", sep="/")
  dir.create(output_dir, recursive=T)

  # Calculates GO enrichment for MCIA
  spls_go_enrichment_results = get_go_enrichment(gene_list = spls_selected_variables_entrezid,
                                                 pvalue_threshold = enrichment_significance_threshold,
                                                 ontology = GO_ontology,
                                                 multiple_test_adjustment = multiple_test_adjustment)
  spls_go_enrichment = spls_go_enrichment_results$data.frame
  spls_TopGOdata = spls_go_enrichment_results$TopGOdata

  if (dim(spls_go_enrichment)[1] > 0) {
    write.table(spls_go_enrichment, file = paste(output_dir, "go_enrichment.txt", sep="/"), sep="\t", row.names = F, quote = F)
    flog.info("GO enrichment on sPLS selected genes returned %d enriched GO terms.", dim(spls_go_enrichment)[1])
    # Plots enrichment for sPLS results
    cluster_and_plot(enrichment_results = spls_go_enrichment,
                     gene_list = spls_results$selected_variables,
                     TopGOdata = spls_TopGOdata,
                     output_dir = output_dir,
                     method = GO_similarity_measure,
                     clustering_method = "pam",
                     ont = GO_ontology,
                     search_universe = GOA_search_universe)
  } else {
    flog.info("GO enrichment on sPLS selected genes gave no results.")
    if (multiple_test_adjustment) {
      flog.info("Try disabling multiple test adjustment.")
    }
  }

}

test <- function()
{

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
go_clustering_evaluation <- function (results_folder, method){

  #base_folder = "results/GO_enrichment"
  base_folder = paste(paste(results_folder, method, sep="/"), "functional_analysis", sep="/")

  # Retrieves results for distance to centroid of each similarity metric
  results_binary = read_clustering_results(paste(base_folder, "binary", sep="/"), "binary")
  dist2centroid = results_binary$dist2centroid

  results_UI = read_clustering_results(paste(base_folder, "UI", sep="/"), "UI")
  dist2centroid = merge(dist2centroid, results_UI$dist2centroid)

  results_BC = read_clustering_results(paste(base_folder, "bray-curtis", sep="/"), "BC")
  dist2centroid = merge(dist2centroid, results_BC$dist2centroid)

  results_cosine = read_clustering_results(paste(base_folder, "cosine", sep="/"), "cosine")
  dist2centroid = merge(dist2centroid, results_cosine$dist2centroid)

  #results_binary_intra = read_clustering_results(paste(method, "binary_gene_list", sep="/"), "binary_intra")
  #dist2centroid = merge(dist2centroid, results_binary_intra$dist2centroid)

  #results_UI_intra = read_clustering_results(paste(method, "UI_gene_list", sep="/"), "UI_intra")
  #dist2centroid = merge(dist2centroid, results_UI_intra$dist2centroid)

  #results_BC_intra = read_clustering_results(paste(method, "bray-curtis_gene_list", sep="/"), "BC_intra")
  #dist2centroid = merge(dist2centroid, results_BC_intra$dist2centroid)

  #results_cosine_intra = read_clustering_results(paste(method, "cosine_gene_list", sep="/"), "cosine_intra")
  #dist2centroid = merge(dist2centroid, results_cosine_intra$dist2centroid)

  results_jiang = read_clustering_results(paste(base_folder, "Jiang", sep="/"), "Jiang")
  dist2centroid = merge(dist2centroid, results_jiang$dist2centroid)

  results_lin = read_clustering_results(paste(base_folder, "Lin", sep="/"), "Lin")
  dist2centroid = merge(dist2centroid, results_lin$dist2centroid)

  results_rel = read_clustering_results(paste(base_folder, "Rel", sep="/"), "Schliker")
  dist2centroid = merge(dist2centroid, results_rel$dist2centroid)

  results_resnik = read_clustering_results(paste(base_folder, "Resnik", sep="/"), "Resnik")
  dist2centroid = merge(dist2centroid, results_resnik$dist2centroid)

  results_wang = read_clustering_results(paste(base_folder, "Wang", sep="/"), "Wang")
  dist2centroid = merge(dist2centroid, results_wang$dist2centroid)


  # Plots boxplot
  dist2centroid_pivot = melt(dist2centroid, id.vars='GO', measure.vars=names(dist2centroid)[2:10])
  png(paste(base_folder,"boxplot_dist2centroid.png", sep="/"))
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
    get_clustering_metrics("Jiang", results_jiang$data),
    get_clustering_metrics("Lin", results_lin$data),
    get_clustering_metrics("Rel", results_rel$data),
    get_clustering_metrics("Resnik", results_resnik$data),
    get_clustering_metrics("Wang", results_wang$data)
  )))


  # Plots barplot with number of clusters
  png(paste(base_folder,"count_clusters.png", sep="/"))
  ggplot(clustering_metrics, aes(x = factor(measure, levels=measure), y = clusters)) +
    geom_bar(stat = "identity", position=position_dodge(width = 0.9), width=0.5) +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5)) +
    ylab("# of clusters") +
    xlab("Measures")
  dev.off()

  clustering_metrics_pivot = melt(clustering_metrics, id.vars='measure', measure.vars=names(clustering_metrics)[4:8])
  clustering_metrics_pivot$value = as.numeric(clustering_metrics_pivot$value)
  clustering_metrics_pivot$measure = factor(clustering_metrics_pivot$measure, levels=clustering_metrics_pivot$measure)

  png(paste(base_folder,"descriptive_stats_clusters.png", sep="/"))
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
  plot_triple_venn(wang_set, resnik_set, cosine_set, c("Wang", "Resnik", "Cosine"), base_folder, "venn_cosine_wang_resnik.png")

  UI_set = unique(results_UI$data$cluster)
  plot_triple_venn(wang_set, resnik_set, UI_set, c("Wang", "Resnik", "UI"), base_folder, "venn_UI_wang_resnik.png")

  plot_triple_venn(wang_set, resnik_set, binary_set, c("Wang", "Resnik", "binary"), base_folder, "venn_binary_wang_resnik.png")

  BC_set = unique(results_BC$data$cluster)
  plot_triple_venn(wang_set, resnik_set, BC_set, c("Wang", "Resnik", "BC"), base_folder, "venn_BC_wang_resnik.png")


}


# Reads cluster results
read_clustering_results <- function(folder, measure){

  #base_folder = "results/GO_enrichment"
  #folder = paste(base_folder, folder, sep="/")
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
