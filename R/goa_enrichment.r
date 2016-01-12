#library(reshape)
#library(biomaRt)
#library(topGO)
#library(GOSemSim)
#library(cluster)


# Gene Ontology Annotations for Human genes and for all Uniprot. No electronically inferred annotations (IEA)
GOA_HUMAN_ANNOTATIONS_URL = "http://geneontology.org/gene-associations/gene_association.goa_human.gz"
GOA_UNIPROT_ANNOTATIONS_URL = "http://geneontology.org/gene-associations/gene_association.goa_uniprot_noiea.gz"
GOA_ANNOTATIONS_FOLDER = get_package_folder("inst/GO/GOA/")

# Reads GOA resource
read_raw_goa <- function(search_universe="human")
{
  dir.create(GOA_ANNOTATIONS_FOLDER, recursive=T)
  if (search_universe == "human"){
    if (! exists("human_goa_raw"))
    {
      file = basename(GOA_HUMAN_ANNOTATIONS_URL)
      file = paste(GOA_ANNOTATIONS_FOLDER, file, sep="")
      download.file(GOA_HUMAN_ANNOTATIONS_URL, file)
      human_goa_raw <<- read.table(file, header = F, comment.char = "!", sep="\t", stringsAsFactors=FALSE, fill = TRUE)
      human_goa_raw <<- human_goa_raw[ , c(3,5)]
      #colnames(human_goa_raw) = c("Gene", "GO")
    }
    colnames(human_goa_raw) = c("Gene", "GO")
    goa_raw <<- human_goa_raw
  } else if (search_universe == "uniprot") {
    if (! exists("uniprot_goa_raw"))
    {
      file = basename(GOA_UNIPROT_ANNOTATIONS_URL)
      file = paste(GOA_ANNOTATIONS_FOLDER, file, sep="")
      download.file(GOA_UNIPROT_ANNOTATIONS_URL, file)
      uniprot_goa_raw <<- read.table(file, header = F, comment.char = "!", sep="\t", stringsAsFactors=FALSE, fill=TRUE)
      uniprot_goa_raw <<- uniprot_goa_raw[ , c(3,5)]
      #colnames(uniprot_goa_raw) = c("Gene", "GO")
    }
    colnames(uniprot_goa_raw) = c("Gene", "GO")
    goa_raw <<- uniprot_goa_raw
  }

  goa_raw
}

# Downloads, parses and formats GOA annotations to match TopGO expected format.
prepare_goa_annotations <- function()
{
  # Reads GOA only for human genes
  goa_raw = read_raw_goa()

  # Pivots table on gene name
  goa = aggregate(data=goa_raw, GO~Gene,c)
  goa$GO = sapply(goa$GO, FUN = function(x) paste(x, collapse=", "))

  # Maps gene identifiers
  map_entrez_hgnc = query_biomart(attributes=c("entrezgene","hgnc_symbol"), filters=c("hgnc_symbol"), values=goa$Gene)

  # Removes HGNC genes not having a map to Entrez
  goa = goa[ goa$Gene %in% map_entrez_hgnc$hgnc_symbol ,]

  # Merges all annotations
  goa = merge(x = goa, y = map_entrez_hgnc, by.x = "Gene", by.y = "hgnc_symbol", all.x = F)

  # Writes to a file the format expected
  goa = goa[ , c("entrezgene", "GO")]
  goa_annotations_dest = paste(GOA_ANNOTATIONS_FOLDER, "goa_parsed.txt", sep="")
  write.table(goa, file = goa_annotations_dest, row.names = F, quote = F, sep = "\t")

  goa_annotations_dest
}

# Returns the semantic similarity score between the two terms
get_semantic_similarity <- function(term1, term2, measure = "Wang", gene_list = NULL, ont="BP", search_universe = NULL)
{
  res <- try(mgoSim(term1, term2, organism="human", measure = measure, ont=ont),silent = TRUE)
  if (class(res) == "try-error" | is.na(res)){
    # Some terms are not found allegedly due to different versions of GO
    res = 0
  }
  res
}

# Returns the functional similarity score between two terms based on the asymmetric binary  or UI distance.
# count(xor) / count(union)
get_functional_similarity <- function(term1, term2, measure = "UI", gene_list = NULL, ont = NULL, search_universe="human")
{
  goa = read_raw_goa(search_universe)
  # if a gene list is provided it reduces the search universe to that list
  if (!is.null(gene_list)){
    goa = goa[goa$Gene %in% gene_list, ]
  }
  term1_genes = unique(goa[goa$GO == term1, c("Gene")])
  term2_genes = unique(goa[goa$GO == term2, c("Gene")])

  # calculates the union
  term_union = union(term1_genes, term2_genes)
  if (length(term_union) == 0 | length(term1_genes) == 0 | length(term2_genes) == 0){
    # When no genes associated to any term they are disimilar
    similarity = 0
  } else {
    # calculates the intersection
    term_intersection = intersect(term1_genes, term2_genes)
    # calculates the different similarity metrics
    if (measure == "binary"){
      term_xor = term_union[!term_union %in% term_intersection]
      similarity = length(term_xor) / length(unique(goa$Gene))
    } else if (measure == "UI"){
      similarity = length(term_intersection) / length(term_union)
    } else if (measure == "bray-curtis"){
      similarity = (2*length(term_intersection)) / (length(term1_genes) + length(term2_genes))
    } else if (measure == "cosine"){
      similarity = length(term_intersection) / (sqrt(length(term1_genes) * length(term2_genes)))
    }
  }
  #term1 = "GO:0000978" term2="GO:0001077"
  similarity
}

# Returns the relative size in percentage of the GO term by the number of genes associatd to it.
get_goa_size <- function(go_term, search_universe="human")
{
  goa_raw = read_raw_goa(search_universe)
  if (! exists("max_goa_size"))
  {
    max_goa_size <<- max(table(unique(goa_raw)$GO))
  }
  (length(unique(goa_raw[goa_raw$GO == go_term , 1])) / max_goa_size)
}

# Calculates the enriched GO terms for a list of genes based on GOA
# Uses the Fisher test provided by TopGO
get_go_enrichment  <- function (gene_list, pvalue_threshold=0.01, ontology="BP", multiple_test_adjustment = FALSE)
{
  flog.info("Computes GO enrichment")

  # Prepares data
  geneID2GO <- readMappings(file = prepare_goa_annotations())
  geneNames <- names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% gene_list))
  names(geneList) <- geneNames
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

  # Performs enrichment test
  test.Fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.Fisher)
  #test.KS <- new("classicScore", testStatistic = GOKSTest, name = "Kolmogorov-Smirnov test")
  #resultKS <- getSigGroups(GOdata, test.KS)

  # The subgraph induced by the top 5 GO terms identified by the elim algorithm for scoring GO terms for
  # enrichment. Rectangles indicate the 5 most significant terms. Rectangle color represents the relative significance,
  # ranging from dark red (most significant) to bright yellow (least significant). For each node, some basic information
  # is displayed. The first two lines show the GO identifier and a trimmed GO name. In the third line the raw p-value
  # is shown. The forth line is showing the number of significant genes and the total number of genes annotated to
  # the respective GO term.
  #showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')

  results <- GenTable(GOdata, classicFisher = resultFisher,
              #classicKS = resultKS, #elimKS = resultKS.elim,
              orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(names(resultFisher@score)))



  if (multiple_test_adjustment){
    # Multiple test correction
    results$pvalue = p.adjust(results$classicFisher, method = "fdr")
  } else {
    results$pvalue = results$classicFisher
  }
  # removes not significant results (keeps those significant)
  results = results[results$pvalue <= pvalue_threshold, ]

  flog.info("Finished computing GO enrichment")

  flog.info("Computes size of GO terms, based on Uniprot-GOA annotations")
  # Calculates frequency (size) of GO terms according to GOA (we use the whole uniprot to calculate size)
  sizes = as.vector(sapply(as.character(results$GO.ID), FUN = get_goa_size))
  flog.info("Finished computing size of GO terms")

  # Returns results
  list(
    TopGOdata = GOdata,
    data.frame = data.frame(
      GO=as.character(results$GO.ID),
      pvalue=results$pvalue,
      name=results$Term,
      annotated_genes=results$Annotated,
      found_genes=results$Significant,
      expected_genes=results$Expected,
      size = sizes,
      row.names=NULL,
      stringsAsFactors = F)
    )
}

# Calculates distances between GO terms based on the semantic similarity
get_go_distance_matrix <- function (enrichment_results, measure = "Wang", gene_list = NULL, ont="BP", search_universe="human")
{
  # Creates a matrix with distances (i.e.: relations) between GO terms considering the genes associated between them
  go_distance_matrix <<- data.frame(matrix(0, nrow = length(enrichment_results$GO), ncol = length(enrichment_results$GO)), stringsAsFactors = F)
  rownames(go_distance_matrix) = enrichment_results$GO
  colnames(go_distance_matrix) = enrichment_results$GO
  if (measure %in% c("UI", "binary", "bray-curtis", "cosine")){
    FUN = get_functional_similarity
  } else {
    FUN = get_semantic_similarity
  }

  # calculates similarity on pairwise comparisons
  pairwise_go_combinations = combn(enrichment_results$GO, 2)
  apply(pairwise_go_combinations, 2, FUN = function(x){
    go_distance_matrix[x[1], x[2]] <<- FUN(term1=x[1], term2=x[2], measure=measure, gene_list=gene_list, ont=ont, search_universe=search_universe)
    })

  # converts similarity matrix into distance matrix
  go_distance_matrix = 1 - go_distance_matrix
  as.dist(t(go_distance_matrix)) # we transpose the matrix as combn calculates similarities on the upper matrix
}

# Hierarchical clustering
hclust_clustering <- function(distance, distance_threshold){
  # Clustering with hclust and a static distance threshold
  clustering = hclust(distance)
  #plot(clustering)
  clusters = cutree(clustering, h=distance_threshold)
  #rect.hclust(fit, h=0.6, border="red")

  clusters
}

# PAM and silhouette clustering
pam_clustering <- function(distance){
  # Clustering with PAM and silhouette
  max_clusters = min(c(30, attr(distance, "Size")-1)) # checks if we have < 30 elements

  flog.info("Number of GO terms %d", attr(distance, "Size")-1)
  flog.info("Looking best number of clusters between %d and %d", 2, max_clusters)
  sil_width = numeric(max_clusters)
  for (k in 2:max_clusters){
    sil_width[k] = pam(distance, k)$silinfo$avg.width
  }
  k_best <- which.max(sil_width[2:max_clusters])
  clustering = pam(distance, k_best)
  clusters = clustering$clustering
  flog.info("Best number of clusters found of %d", k_best)
  clusters
}

# MDS
multidimensional_scaling <- function(distance){
  number_elements = attr(distance, "Size")
  flog.info("Applying MDS to %d elements", number_elements)
  if (number_elements < 2){
    x = rep(0, number_elements)
    y = rep(0, number_elements)
  } else if (number_elements == 2){
    x = c(-distance[1]/2, distance[1]/2)
    y = rep(0, number_elements)
  } else {
    fit <- cmdscale(distance, eig = TRUE, k = 2)
    x <- fit$points[, 1]
    if (dim(fit$points)[2] > 1){
      y = fit$points[, 2]
    } else {
      y = rep(0, length(x))
      names(y) = names(x)
    }
  }

  list(x=x, y=y)
}

# Plots GO terms according to the enrichment results
cluster_and_plot <- function (enrichment_results, gene_list, TopGOdata, method="Wang", search_universe="human", output_dir=".", distance_threshold = 0.7, frequency_threshold=0.05, clustering_method = "hclust", ont="BP")
{

  # Calculates the dissimilarity between GO terms
  if (method %in% c("Resnik", "Lin", "Rel", "Jiang", "Wang", "UI", "binary", "bray-curtis", "cosine", "all") ) {

    # Prepares table
    enrichment_results = cbind(enrichment_results, cluster=rep("", dim(enrichment_results)[1]), cluster_name=rep("", dim(enrichment_results)[1]))
    enrichment_results$cluster = as.character(enrichment_results$cluster)
    enrichment_results$cluster_name = as.character(enrichment_results$cluster_name)

    # Removes all terms with a frequency in GOA >5%
    #enrichment_results = enrichment_results[enrichment_results$size <= frequency_threshold, ]

    if (method == "all"){
      methods = c("Resnik", "Lin", "Rel", "Jiang", "Wang", "UI", "binary", "bray-curtis", "cosine")
    } else {
      methods = c(method)
    }

    parent_output_dir = output_dir
    for (methodN in methods){

      output_dir = paste(parent_output_dir, methodN, sep="/")
      dir.create(output_dir)

      # Obtains the distance matrix between terms
      flog.info("Computing pairwise distance matrix on enriched GO terms with method %s", methodN)
      if (search_universe %in% c("human", "uniprot")){
        go_distance = get_go_distance_matrix(enrichment_results, measure = methodN, ont=ont, search_universe=search_universe)
      } else if (search_universe == "gene_list") {
        # similarity calculations are based on GOA associations limited to those genes in the gene list
        # only applicable to binary, UI, Bray-Curtis and cosine distances
        go_distance = get_go_distance_matrix(enrichment_results, measure = methodN, gene_list = gene_list, ont=ont, search_universe=search_universe)
      }
      flog.info("Finished computing distances")

      # Clustering
      flog.info("Computing clustering on enriched GO terms with method %s", clustering_method)
      if (clustering_method == "hclust"){
        clusters = hclust_clustering(go_distance, distance_threshold)
      } else if (clustering_method == "pam"){
        clusters = pam_clustering(go_distance)
      } else {
        stop ("Clustering method not supported")
      }
      flog.info("Finished computing clustering")

      # Select cluster representatives
      flog.info("Selection of cluster representatives")
      for (cluster in unique(sort(clusters))){
        group = names(clusters)[clusters %in% cluster]
        # Selects as the cluster representative the most significantly associated term with a frequency <= 5% if any
        candidates = enrichment_results[enrichment_results$GO %in% group & enrichment_results$size <= frequency_threshold, c("pvalue")]
        if (length(candidates)==0){
          representative_idx = which.min(enrichment_results[enrichment_results$GO %in% group, c("pvalue")])
        } else {
          representative_idx = which.min(enrichment_results[enrichment_results$GO %in% group & enrichment_results$size <= frequency_threshold, c("pvalue")])
        }
        enrichment_results[enrichment_results$GO %in% group, c("cluster")] = group[representative_idx]
        enrichment_results[enrichment_results$GO %in% group, c("cluster_name")] = enrichment_results[enrichment_results$GO == group[representative_idx], c("name")]
      }
      cluster_representatives = enrichment_results[enrichment_results$GO == enrichment_results$cluster, ]

      # Multi Dimensional Scaling on all terms
      flog.info("Computing Multi Dimensional Scaling")
      mds_results = multidimensional_scaling(go_distance)
      enrichment_results = cbind(enrichment_results, mds_results$x)
      enrichment_results = cbind(enrichment_results, mds_results$y)

      # Stores the results of clustering and MDS
      write.table(enrichment_results, file=paste(output_dir, "enrichment_results.txt", sep="/"), sep="\t", quote = F, row.names = F)

      # Stores Euclidean distance to centroid to evaluate dispersion
      avg_x = mean(mds_results$x)
      avg_y = mean(mds_results$y)
      distances_to_centroid = sqrt((mds_results$x-avg_x)^2 + (mds_results$y-avg_y)^2)
      write.table(data.frame(GO=names(distances_to_centroid), distance_to_centroid=distances_to_centroid), file=paste(output_dir, "distances_to_centroid.txt", sep="/"), sep="\t", quote = F, row.names = F)

      # Multi Dimensional Scaling just on clusters
      go_distance_m = as.matrix(go_distance)
      cluster_distance = as.dist(go_distance_m[rownames(go_distance_m) %in% cluster_representatives$GO, colnames(go_distance_m) %in% cluster_representatives$GO])
      mds_results = multidimensional_scaling(cluster_distance)
      x = mds_results$x
      y = mds_results$y
      cluster_representatives = cbind(cluster_representatives, x)
      cluster_representatives = cbind(cluster_representatives, y)

      flog.info("Finished MDS")

      # Stores the results of clustering and MDS on representatives
      write.table(cluster_representatives, file=paste(output_dir, "cluster_representatives.txt", sep="/"), sep="\t", quote = F, row.names = F)

      # Plots
      flog.info("Plots functional analysis results")
      plot_scatter(cluster_representatives, output_dir)
      plot_table(cluster_representatives = cluster_representatives, output_dir)
      plot_treemap(enrichment_results = enrichment_results, output_dir)
      plot_graph(cluster_representatives = cluster_representatives, GOdata = TopGOdata, output_dir)
      flog.info("Finished plotting")
    }
  } else {
    # by default runs Wang semantic similarity
    stop("Non supported method.")
  }
}
