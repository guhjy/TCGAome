library(ReactomePA)
library(RTCGAToolbox)
library(mixOmics)
library(omicade4)


# COlors to use in plots
red500 = rgb(as.integer(as.hexmode("f4")), as.integer(as.hexmode("43")), as.integer(as.hexmode("36")), maxColorValue=255)
green500 = rgb(as.integer(as.hexmode("4c")), as.integer(as.hexmode("af")), as.integer(as.hexmode("50")), maxColorValue=255)
indigo500 = rgb(as.integer(as.hexmode("3f")), as.integer(as.hexmode("51")), as.integer(as.hexmode("b5")), maxColorValue=255)
amber500 = rgb(as.integer(as.hexmode("ff")), as.integer(as.hexmode("c1")), as.integer(as.hexmode("07")), maxColorValue=255)
yellow500 = rgb(as.integer(as.hexmode("ff")), as.integer(as.hexmode("eb")), as.integer(as.hexmode("3b")), maxColorValue=255)
purple500 = rgb(as.integer(as.hexmode("67")), as.integer(as.hexmode("3a")), as.integer(as.hexmode("b7")), maxColorValue=255)
grey500 = rgb(as.integer(as.hexmode("9e")), as.integer(as.hexmode("9e")), as.integer(as.hexmode("b9")), maxColorValue=255)
brown500 = rgb(as.integer(as.hexmode("79")), as.integer(as.hexmode("55")), as.integer(as.hexmode("48")), maxColorValue=255)
cyan500 = rgb(as.integer(as.hexmode("00")), as.integer(as.hexmode("bc")), as.integer(as.hexmode("d4")), maxColorValue=255)
lime500 = rgb(as.integer(as.hexmode("cd")), as.integer(as.hexmode("dc")), as.integer(as.hexmode("39")), maxColorValue=255)
black = rgb(as.integer(as.hexmode("00")), as.integer(as.hexmode("00")), as.integer(as.hexmode("00")), maxColorValue=255)




###
## Dataset descriptive analysis
###
descriptive.analysis <- function(X, Y, Z){
  output_dir = paste(RESULTS_FOLDER, "descriptive_analysis", sep="/")
  dir.create(output_dir)
  
  plot_samples_barplot(Z, output_dir=output_dir, file_name="samples.barplot.png")
}


###
## Principal Component Analysis
###
pca.analysis <- function(X, Y, Z){
  
  output_dir = paste(RESULTS_FOLDER, "PCA", sep="/")
  dir.create(output_dir)
  
  # PCA on the protein expression
  Y.pca <- prcomp(Y, center = TRUE, scale. = TRUE)
  plot_simple_scatter(matrix=Y.pca, groups=Z$Tumor_type, output_dir=output_dir, file_name="PCA.Y.png")
  
  # PCA on the gene expression
  X.pca <- prcomp(X,
                  center = TRUE,
                  scale. = TRUE)
  plot_simple_scatter(matrix=X.pca, groups=Z$Tumor_type, output_dir=output_dir, file_name="PCA.X.png")
  
  # Joint PCA
  colnames(X) = paste(colnames(X), ".X")
  colnames(Y) = paste(colnames(Y), ".Y")
  XY.pca <- prcomp(cbind2(X,Y),
                   center = TRUE,
                   scale. = TRUE) 
  plot_simple_scatter(matrix=XY.pca, groups=Z$Tumor_type, output_dir=output_dir, file_name="PCA.XY.png")
}

###
## Hierarchical clustering
###
hclust.analysis <- function(X, Y, Z){
  
  output_dir = paste(RESULTS_FOLDER, "hclust", sep="/")
  dir.create(output_dir)
  
  col.tumor.type = Z$Tumor_type
  col.tumor.type[col.tumor.type == "BRCA"] <- red500
  col.tumor.type[col.tumor.type == "OV"] <- green500
  
  colnames(X) = paste(colnames(X), ".X")
  colnames(Y) = paste(colnames(Y), ".Y")
  
  # Plots heatmaps
  plot_heatmap(matrix=as.matrix(X), colors=col.tumor.type, output_dir=output_dir, file_name="hclust.X.png")
  plot_heatmap(matrix=as.matrix(Y), colors=col.tumor.type, output_dir=output_dir, file_name="hclust.Y.png")
  plot_heatmap(matrix=as.matrix(cbind2(X,Y)), colors=col.tumor.type, output_dir=output_dir, file_name="hclust.XY.png")
}

###
### Regularized Canonical COrrelation Analysis (rCCA)
###
rcaa.analysis <- function(X, Y, Z)
{
  
  output_dir = paste(RESULTS_FOLDER, "rCCA", sep="/")
  dir.create(output_dir)
  
  # Calculates regularization factors with CV (M-folds and leave-one-out methods raise an error)
  # Error in chol.default(Cxx) : 
  # the leading minor of order 407 is not positive definite
  #grid1 <- seq(0, 0.2, length = 51)
  #grid2 <- seq(0.0001, 0.2, length = 51)
  cv_score <- tune.rcc(X, Y)  #, grid1 = grid1, grid2 = grid2, plt = FALSE) #, validation = 'loo')
  # Computes rCCA
  rcca_results <- rcc(X, Y, lambda1=cv_score$opt.lambda1, lambda2=cv_score$opt.lambda2)
  
  # plots the individuals
  col.tumor.type = Z$Tumor_type
  col.tumor.type[col.tumor.type == "BRCA"] <- red500
  col.tumor.type[col.tumor.type == "OV"] <- green500
  mixomics_plot_individuals(data=rcca_results, names=Z$race, colors=col.tumor.type, output_dir=output_dir, file_name="samples.png")
  
  # plots the variables
  mixomics_plot_variables(rcca_results, output_dir=output_dir, file_name="variables.png")
  
  # network of relations between variables
  #network(rcca_results, comp = 1:2, threshold = 0.2)
  
  # heatmap
  #cim(rcca_results, comp = 1:2, xlab = "proteins", ylab = "genes", margins = c(5,6))
  
  # return
  rcca_results
}

###
## Regularized Generalized CCA
###
rgcca.analysis <- function (X, Y, Z)
{
  library(RGCCA)
  
  rgcca.results = rgcca(list(X, Y), tau = 'optimal', ncomp = 3)
  
  # return
  rgcca.results
}


###
### Sparse Partial Least Squares (sPLS)
###
spls.analysis <- function(X, Y, Z, topN=10)
{
  
  output_dir = paste(RESULTS_FOLDER, "mcia", sep="/")
  dir.create(output_dir)
  
  # adds a suffix to gene names to differentiate the dataset of origin, we'll have to remove for posterior analysis
  colnames(X) = paste(colnames(X), ".X", sep="")
  colnames(Y) = paste(colnames(Y), ".Y", sep="")
  
  #spls_result <- spls(as.matrix(X), as.matrix(Y), ncomp = 3, mode = 'regression', keepX=c(50, 50, 50), keepY=c(50, 50, 50))
  spls_result <- spls(as.matrix(X), as.matrix(Y), ncomp = 3, mode = 'regression')
  
  col.tumor.type = Z$Tumor_type
  col.tumor.type[col.tumor.type == "BRCA"] <- red500
  col.tumor.type[col.tumor.type == "OV"] <- green500
  
  # plots the individuals
  mixomics_plot_individuals(data=spls_result, names=Z$Tumor_type, colors=col.tumor.type, output_dir=output_dir, file_name="samples.png")
  # plots the variables
  mixomics_plot_variables(spls_result, output_dir=output_dir, file_name="variables.png")
  # plots variables in a network
  mixomics_plot_network(spls_result, output_dir=output_dir, file_name="network.png")
  # plots heatmap
  mixomics_plot_heatmap(spls_result, output_dir=output_dir, file_name="heatmap.png")
  
  # Select variables
  selected_variables = mixOmics::selectVar(spls_result, comp=3)
  rownames(selected_variables$value.X) = gsub("\\.X", "", rownames(selected_variables$value.X))
  rownames(selected_variables$value.Y) = gsub("\\.Y", "", rownames(selected_variables$value.Y))
  selected_variables$value.X$name.X = rownames(selected_variables$value.X)
  selected_variables$value.Y$name.Y = rownames(selected_variables$value.Y)
  
  
  positive.X = selected_variables$value.X[selected_variables$value.X$value.var.X >=0, ]
  b = dim(positive.X)[1]
  a = b - topN + 1
  positive.X.top = positive.X[a:b , ]
  negative.X = selected_variables$value.X[selected_variables$value.X$value.var.X <0, ]
  b = dim(negative.X)[1]
  a = b - topN + 1
  negative.X.top = negative.X[a:b , ]
  positive.Y = selected_variables$value.Y[selected_variables$value.Y$value.var.Y >=0, ]
  b = dim(positive.Y)[1]
  a = b - topN + 1
  positive.Y.top = positive.Y[a:b , ]
  negative.Y = selected_variables$value.Y[selected_variables$value.Y$value.var.Y <0, ]
  b = dim(negative.Y)[1]
  a = b - topN + 1
  negative.Y.top = negative.Y[a:b , ]
  
  selected_variables_names = c(positive.X.top$name.X, negative.X.top$name.X, positive.Y.top$name.Y, negative.Y.top$name.Y)
  selected_genes = c(positive.X.top$name.X, negative.X.top$name.X)
  selected_proteins = c(positive.Y.top$name.Y, negative.Y.top$name.Y)
  
  # return
  list(results=spls_result, selected_variables=selected_variables_names, selected_genes=selected_genes, selected_proteins=selected_proteins)
}

###
### Co-Inertia Analysis
###
mcia.analysis <- function(X, Y, Z, topN=5, cia.nf=5)
{
  
  output_dir = paste(RESULTS_FOLDER, "mcia", sep="/")
  dir.create(output_dir)
  
  col.histological.type = as.numeric(Z$histological_type)
  col.histological.type[col.histological.type == 1] <- red500
  col.histological.type[col.histological.type == 2] <- green500
  col.histological.type[col.histological.type == 3] <- indigo500
  col.histological.type[col.histological.type == 4] <- amber500
  col.histological.type[col.histological.type == 5] <- yellow500
  col.histological.type[col.histological.type == 6] <- purple500
  
  col.disease.stage = as.numeric(Z$neoplasm_diseasestage)
  col.disease.stage[col.disease.stage == 1] <- red500
  col.disease.stage[col.disease.stage == 2] <- green500
  col.disease.stage[col.disease.stage == 3] <- indigo500
  col.disease.stage[col.disease.stage == 4] <- amber500
  col.disease.stage[col.disease.stage == 5] <- yellow500
  col.disease.stage[col.disease.stage == 6] <- purple500
  col.disease.stage[col.disease.stage == 7] <- brown500
  col.disease.stage[col.disease.stage == 8] <- cyan500
  col.disease.stage[col.disease.stage == 9] <- lime500
  col.disease.stage[col.disease.stage == 10] <- black
  col.disease.stage[is.na(col.disease.stage)] <- grey500
  
  
  #matrices_list = list(t(X.tsg.subset), t(Y))
  #matrices_list = list(t(X), t(Y), t(Z))
  matrices_list = list(as.matrix(t(X)), as.matrix(t(Y)))
  sapply(matrices_list, dim)
  all(apply((x <- sapply(matrices_list, colnames)), 2, function(y)identical(y, x[,1])))
  
  # Multiple Co-Inertia Analysis
  mcia_result <- mcia(matrices_list, nsc=F, cia.nf = cia.nf)
  
  # PLot all results
  mcia_plot(mcia_result, phenotype=Z$Tumor_type, output_dir=output_dir, file_name="visualizations.png")
  
  # Selects top 5 positive and negative associations and plots them, that is 20 variables
  mcia_selected_variables = topVar(mcia_result, topN = topN)
  mcia_selected_genes = gsub("\\.df1", "", as.character(unlist(mcia_selected_variables[, 1:2])))
  mcia_selected_proteins = gsub("\\.df2", "", as.character(unlist(mcia_selected_variables[, 3:4])))
  write.table(mcia_selected_variables, paste(output_dir, "top_variables.txt", sep="/"), quote = FALSE, sep = "\t", row.names = FALSE)
  mcia_selected_variables = as.character(unlist(mcia_selected_variables))
  
  # Plots selected variables
  mcia_plot_variables(mcia_result, mcia_selected_variables, output_dir=output_dir, file_name="topN.variables.png")
  
  # Removes the suffix that identifies a variable as gene or protein
  mcia_selected_variables = gsub("\\.df.", "", mcia_selected_variables)
  
  list(results=mcia_result, selected_variables=mcia_selected_variables, selected_genes=mcia_selected_genes, selected_proteins=mcia_selected_proteins)
}


###
### Results evaluation
###
results.evaluation <- function(spls.results, mcia.results)
{
  
  # TSGs and oncogenes from publication Vogelstein et al. (2013). Cancer genome landscapes. Science
  oncogenes_and_tsg = read.csv("../data/oncogenes_and_tsg.processed.txt", sep="\t")
  oncogenes = as.character(oncogenes_and_tsg[oncogenes_and_tsg$Classification == "Oncogene", 1])
  tsgs = as.character(oncogenes_and_tsg[oncogenes_and_tsg$Classification == "TSG", 1])
  #oncogenes = intersect(colnames(X), oncogenes)
  #tsgs = intersect(colnames(X), tsgs)
  as.character(oncogenes_and_tsg$Gene.Symbol)
  
  
  # Plots Venn diagram
  plot_triple_venn(
    mcia.results$selected_variables,
    spls.results$selected_variables,
    as.character(oncogenes_and_tsg$Gene.Symbol),
    c("MCIA", "sPLS", "Oncogenes and TSGs",),
    output_dir="../results",
    file_name="venn.png"
  )
  
  
  # TSGs in BRCA tumors selected by differential gene expression pan-cancer analysis
  brca_tsg = read.table("../data/BRCA_down_regulated_TSgenes.txt", col.names = c("gene_id", "gene_name"))
  brca_tsg = intersect(colnames(X), brca_tsg$gene_name)
  
  # BRCA TSGs selected by MCIA in the top 5
  known_brca_tsg_selected_by_mcia = intersect(brca_tsg, mcia_selected_variables)
  brca_tsg_variables = paste(brca_tsg, ".df1", sep="")
  # Plots BRCA TSGs against MCIA variable space
  plotVar(mcia_result, brca_tsg_variables, var.col=red500, bg.var.col="grey")
  
  # Oncogenes and TSGs selected by MCIA in the top 5
  known_tsgs_selected_by_mcia = intersect(tsgs, mcia_selected_variables)
  known_oncogenes_selected_by_mcia = intersect(oncogenes, mcia_selected_variables)
  
  # Plots TSGs against MCIA variable space
  plotVar(mcia_result, paste(oncogenes, ".df1", sep=""), var.col=red500, bg.var.col="grey")
  plotVar(mcia_result, paste(tsgs, ".df1", sep=""), var.col=red500, bg.var.col="grey")
  
  
  # Calculate the probability of finding a TSG or oncogene in the top 10 by chance
  trials = 10
  n = dim(X)[2]
  n1 = length(brca_tsg)
  n2 = length(oncogenes)
  
  # p1 = 0.2672634; p2 = 0.02786402
  binomial_coefficient <- function(n, k) exp(lfactorial(n) - lfactorial(k) - lfactorial(n-k))
  p1 = 1 - exp(log(binomial_coefficient(n-n1, trials)) - log(binomial_coefficient(n, trials)))
  p2 = 1 - exp(log(binomial_coefficient(n-n2, trials)) - log(binomial_coefficient(n, trials)))
}


###
### Functional analysis
###
functional.analysis <- function (spls.results, mcia.results)
{
  
  # Loads ensembl biomart for Homo sapiens
  
  
  oncogenes_and_tsg = read.csv("./data/oncogenes_and_tsg.processed.txt", sep="\t")
  
  # get entrezid gene list from genes (it loses one gene on conversion, I guess it is "C11orf75")
  mcia_selected_variables_entrezid = query_biomart(attributes=c("entrezgene"), filters=c("hgnc_symbol"), values=mcia.results$selected_variables)[, 1]
  spls_selected_variables_entrezid = query_biomart(attributes=c("entrezgene"), filters=c("hgnc_symbol"), values=spls.results$selected_variables)[, 1]
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
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/Wang", 
                   method = "Wang",
                   clustering_method = "pam",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/binary", 
                   method = "binary",
                   clustering_method = "pam",
                   search_universe = "human",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/binary_gene_list", 
                   method = "binary",
                   clustering_method = "pam",
                   search_universe = "gene_list",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/binary_uniprot", 
                   method = "binary",
                   clustering_method = "pam",
                   search_universe = "uniprot",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/UI", 
                   method = "UI",
                   clustering_method = "pam",
                   search_universe = "human",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/UI_gene_list", 
                   method = "UI",
                   clustering_method = "pam",
                   search_universe = "gene_list",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/UI_uniprot", 
                   method = "UI",
                   clustering_method = "pam",
                   search_universe = "uniprot",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/bray-curtis", 
                   method = "bray-curtis",
                   clustering_method = "pam",
                   search_universe = "human",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/bray-curtis_gene_list", 
                   method = "bray-curtis",
                   clustering_method = "pam",
                   search_universe = "gene_list",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/bray-curtis_uniprot", 
                   method = "bray-curtis",
                   clustering_method = "pam",
                   search_universe = "uniprot",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/cosine", 
                   method = "cosine",
                   clustering_method = "pam",
                   search_universe = "human",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/cosine_gene_list", 
                   method = "cosine",
                   clustering_method = "pam",
                   search_universe = "gene_list",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/cosine_uniprot", 
                   method = "cosine",
                   clustering_method = "pam",
                   search_universe = "uniprot",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/Resnik", 
                   method = "Resnik",
                   clustering_method = "pam",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/Lin", 
                   method = "Lin",
                   clustering_method = "pam",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
                   TopGOdata = mcia_TopGOdata, 
                   output_dir="results/GO_enrichment/MCIA/Rel", 
                   method = "Rel",
                   clustering_method = "pam",
                   ont=ONTOLOGY)
  
  cluster_and_plot(enrichment_results = mcia_go_enrichment, 
                   gene_list = mcia.results$selected_variables,
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
  library(dplyr)
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
  
  
  library(VennDiagram)
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



###
## Main
###
run.TCGAome <- function(tumor_types){
  
  # Creates and sets the results folder for the current run
  RESULTS_FOLDER <<- get_results_folder()
  
  # Downloads data for RNAseq and RPPA
  matrices = get.data(c("BRCA", "OV"))
  
  # Preprocessing
  preprocessed.matrices = preprocess.data(X=matrices$X, Y=matrices$Y, correlation.thr = 0.7)
  
  # Pre-analysis
  descriptive.analysis(X=preprocessed.matrices$X, Y=preprocessed.matrices$Y, Z=matrices$Z)
  pca.analysis(X=preprocessed.matrices$X, Y=preprocessed.matrices$Y, Z=matrices$Z)
  hclust.analysis(X=preprocessed.matrices$X, Y=preprocessed.matrices$Y, Z=matrices$Z)
  
  # Runs MCIA
  mcia.results = mcia.analysis(X=preprocessed.matrices$X, Y=preprocessed.matrices$Y, Z=matrices$Z, topN=10, cia.nf=5)
  
  # Runs sPLS
  #spls.preprocessed.matrices = preprocess.data(matrices$X, matrices$Y, correlation.thr = 0.5)
  #spls.results = spls.analysis(X=spls.preprocessed.matrices$X, Y=spls.preprocessed.matrices$Y, Z=matrices$Z)
  spls.results = spls.analysis(X=preprocessed.matrices$X, Y=preprocessed.matrices$Y, Z=matrices$Z)
  
  # Runs RGCCA
  rgcca.results = rgcca.analysis(X=preprocessed.matrices$X, Y=preprocessed.matrices$Y, Z=matrices$Z)
  
  # Runs rCCA
  #preprocessed.matrices = preprocess.data(matrices$X, matrices$Y, correlation.thr = 0.5)
  rcaa.results = rcaa.analysis(preprocessed.matrices$X, preprocessed.matrices$Y, Z=matrices$Z)
  
  # Evaluates results
  results.evaluation(mcia.results = mcia.results, spls.results = spls.results)
  
  # Performs pathway enrichment analysis
  functional.analysis(mcia.results = mcia.results, spls.results = spls.results)
  
  
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



