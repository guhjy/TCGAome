
# Gets the folder to store GOA resources
get_goa_folder <- function() {
    get_package_folder("inst/GO/GOA/")
}

# Gets the folder to store GOA resources
get_hpo_folder <- function() {
    get_package_folder("inst/HPO/")
}


#' Reads Gene Ontology Annotations (GOA) for Human genes and for all Uniprot. No electronically inferred annotations (IEA)
#' @param search_universe Two possible values "human" or "uniprot". The first provides associations to GO of human genes,
#' while the second considers all genes registered in Uniprot. [default: "human"]
#' @param ontology The ontology within Gene Ontology: "BP", "CC" or "MF". [default: "BP"]
#' @param goa_human_annotations_url The URL to human GOA. [default: "http://geneontology.org/gene-associations/gene_association.goa_human.gz"]
#' @param goa_uniprot_annotations_url The URL to Uniprot GOA. [default: "http://geneontology.org/gene-associations/gene_association.goa_uniprot_noiea.gz"]
#'
#' @keywords TCGAome
#' @export
#' @examples
#' goa = load_goa()
load_goa <- function(search_universe = "human",
                     ontology = "BP",
                     goa_uniprot_annotations_url = "http://geneontology.org/gene-associations/gene_association.goa_uniprot_noiea.gz") {

    ## get and create working folder
    goa_annotations_folder <- TCGAome::get_goa_folder()
    goa = NULL

    if (search_universe == "human") {
        ## Loads stored attribute if any
        goa <- attr(load_goa, "human_goa")
        stored_ontology <- attr(load_goa, "ontology")
        if (is.null(goa) | ontology != stored_ontology) {
            ## Downloads and stores as attribute human GOA
            futile.logger::flog.debug("Loading human GOA for the first time")
            uniKeys <- keys(org.Hs.eg.db, keytype="SYMBOL")
            cols <- c("GOALL")
            goa_raw <- select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="SYMBOL")
            goa_raw <- goa_raw[goa_raw$ONTOLOGYALL == ontology, ]
            goa_raw <- goa_raw[, c(1, 2)]
            colnames(goa_raw) <- c("Gene", "Term")
            goa <- new("GeneAnnotations", raw_annotations = goa_raw, name="GOA-Human")
            attr(load_goa, "human_goa") <<- goa
            attr(load_goa, "ontology") <<- ontology
        } else {
            futile.logger::flog.debug("Loading cached human GOA")
        }
    } else if (search_universe == "uniprot") {
        ## Loads stored attribute if any
        goa <- attr(load_goa, "uniprot_goa")
        if (is.null(goa)) {
            ## Downloads and stores as attribute uniprot GOA
            futile.logger::flog.debug("Loading uniprot GOA for the first time")
            file <- basename(goa_uniprot_annotations_url)
            file <- paste(goa_annotations_folder, file, sep = "")
            download.file(goa_uniprot_annotations_url, file, quiet = TRUE)
            goa_raw <- read.table(file, header = F, comment.char = "!", sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
            goa_raw <- goa_raw[, c(3, 5)]
            colnames(goa_raw) <- c("Gene", "Term")
            ## Preprocess the GOA table, maps gene ids to entrez gene ids
            goa <- new("GeneAnnotations", raw_annotations = goa_raw, name="GOA-Uniprot")
            attr(load_goa, "uniprot_goa") <<- goa
        } else {
            futile.logger::flog.debug("Loading cached uniprot GOA")
        }
    } else {
        futile.logger::flog.error("Non supported search universe [%s]", search_universe)
    }

    return(goa)
}

#' Reads KEGG pathways for Human genes and for all Uniprot. No electronically inferred annotations (IEA)
#'
#' @keywords TCGAome
#' @export
#' @examples
#' goa = load_kegg()
load_kegg <- function() {

    kegg = NULL

    ## Loads stored attribute if any
    kegg <- attr(load_kegg, "human_kegg")
    if (is.null(kegg)) {
        ## Downloads and stores as attribute human GOA
        futile.logger::flog.debug("Loading human KEGG for the first time")
        uniKeys <- keys(org.Hs.eg.db, keytype="SYMBOL")
        cols <- c("PATH")
        kegg_raw <- select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="SYMBOL")
        kegg_raw <- kegg_raw[, c(1, 2)]
        colnames(kegg_raw) <- c("Gene", "Term")
        kegg <- new("GeneAnnotations", raw_annotations = kegg_raw, name="KEGG-Human")
        attr(load_kegg, "human_kegg") <<- kegg
    } else {
        futile.logger::flog.debug("Loading cached human KEGG")
    }

    return(kegg)
}

#' Reads OMIM diseases for Human genes.
#'
#' @keywords TCGAome
#' @export
#' @examples
#' goa = load_omim()
load_omim <- function() {

    omim = NULL

    ## Loads stored attribute if any
    omim <- attr(load_omim, "omim")
    if (is.null(omim)) {
        ## Downloads and stores as attribute human OMIM
        futile.logger::flog.debug("Loading human OMIM for the first time")
        uniKeys <- keys(org.Hs.eg.db, keytype="SYMBOL")
        cols <- c("OMIM")
        omim_raw <- select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="SYMBOL")
        omim_raw <- omim_raw[, c(1, 2)]
        colnames(omim_raw) <- c("Gene", "Term")
        omim <- new("GeneAnnotations", raw_annotations = omim_raw, name="OMIM-Human")
        attr(load_omim, "human_omim") <<- omim
    } else {
        futile.logger::flog.debug("Loading cached human OMIM")
    }

    return(omim)
}

#' Reads Human Phenotype Ontology (HPO) annotations.
#' @param hpo_annotations_url The URL to HPO annotations to genes. [default: "http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_TYPICAL_FEATURES_phenotype_to_genes.txt"]
#'
#' @keywords TCGAome
#' @export
#' @examples
#' hpo = load_hpo()
load_hpo <- function(
    hpo_annotations_url =
        "http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_TYPICAL_FEATURES_phenotype_to_genes.txt") {

    ## get and create working folder
    hpo_annotations_folder <- TCGAome::get_hpo_folder()

    ## Loads stored attribute if any
    hpo <- attr(load_hpo, "hpo")
    if (is.null(hpo)) {
        ## Downloads and stores as attribute human GOA
        futile.logger::flog.debug("Loading HPO for the first time")
        file <- basename(hpo_annotations_url)
        file <- paste(hpo_annotations_folder, file, sep = "")
        download.file(hpo_annotations_url, file, quiet = TRUE)
        hpo_raw <- read.table(file, header = F, comment.char = "#", sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
        hpo_raw <- hpo_raw[, c(4, 1)]
        colnames(hpo_raw) <- c("Gene", "Term")
        ## Preprocess the HPO table, maps gene ids to entrez gene ids
        hpo <- new("GeneAnnotations", raw_annotations = hpo_raw, name="HPO")
        attr(load_goa, "hpo") <<- hpo
    } else {
        futile.logger::flog.debug("Loading cached HPO")
    }

    return(hpo)
}


#' Computes the enrichment of a gene list into a given annotation resource (GOA and HPO are supported so far).
#' The enrichment is computed with the Fisher's exact test, multiple test
#' correction is optional.
#' @param gene_list The list of genes (entrez ids) for which we want to calculate the enrichment
#' @param annotation One of ["goa", "hpo"]
#' @param pvalue_threshold The threshold that determines significant results [default: 0.01]
#' @param adjustment_method Any method supported by p.adjust() [default: none]
#'
#' @keywords TCGAome
#' @export
#' @examples
#' size = get_enrichment(c("57763", "81611", "367", "84002", "578", "51161", "840", "891",
#' "595", "8318", "1001", "1291", "1408", "1857", "30846", "2099",
#' "2146", "26234", "2309", "8833", "2931", "3192", "3303", "8821",
#' "9929", "3791", "122773", "3930", "254531", "9833", "10403",
#' "4790", "5163", "5241", "644591", "64284", "65059", "5925",
#' "6584", "4089", "6774", "3925", "100500808", "6850", "7249",
#' "123016", "7874", "7515", "10413", "4904", "7542", "27332",
#' "57116", "116412"), "goa")
get_enrichment <- function(gene_list, annotation, pvalue_threshold = 0.01, adjustment_method = "none") {

    futile.logger::flog.info("Computing enrichment with Fisher's test...")

    ## Prepares data
    if (annotation == "goa") {
        resource <- TCGAome::load_goa()
    } else if (annotation == "hpo") {
        resource <- TCGAome::load_hpo()
    } else {
        futile.logger::flog.error("Non supported annotation [%s]", annotation)
        return()
    }

    all_genes <- resource@gene2term$Gene
    selected_genes <- gene_list
    all_terms <- resource@term2gene$Term

    ## Computes the Fisher test for every GO term
    results = sapply(all_terms, function(X) {
        associated_genes = as.character(unlist(resource@term2gene[resource@term2gene$Term == X, ]$Gene))
        fisher.test(matrix(c(
            length(intersect(selected_genes, associated_genes)),
            length(associated_genes[! associated_genes %in% selected_genes]),
            length(selected_genes[! selected_genes %in% associated_genes]),
            length(all_genes[! all_genes %in% union(associated_genes, selected_genes)])
        ), nrow = 2, ncol = 2),
        alternative = "greater")$p.value
    })

    # Multiple test correction
    futile.logger::flog.info("Applying [%s] multiple test correction method", adjustment_method)
    results_adjusted <- p.adjust(results, method = adjustment_method)

    # removes not significant results (keeps those significant)
    results_significant <- results_adjusted[results_adjusted <= pvalue_threshold]
    futile.logger::flog.info("Found %s significantly enriched terms", length(results_significant))

    futile.logger::flog.info("Computing relative frequency of terms in the annotations")
    # Calculates frequency (size) of GO terms according to GOA (we use the whole uniprot to calculate size)
    sizes <- as.double(sapply(names(results_significant), FUN = function(term){
        TCGAome::get_term_size(resource, term)
        }))

    # Returns results
    return(data.frame(
        Term = names(results_significant),
        pvalue = as.double(results_significant),
        size = sizes,
        row.names = NULL, stringsAsFactors = F))

        #name = results$Term,
        #annotated_genes = results$Annotated,
        #found_genes = results$Significant,
        #expected_genes = results$Expected,
}



# Calculates the enriched GO terms for a list of genes based on GOA Uses the Fisher test provided by TopGO
get_go_enrichment <- function(gene_list, pvalue_threshold = 0.01, ontology = "BP", multiple_test_adjustment = FALSE) {
    futile.logger::flog.info("Computes GO enrichment")

    # Prepares data
    goa = TCGAome::load_goa()
    goa <- goa[, c("entrezgene", "GO")]
    goa_annotations_dest <- paste(TCGAome::get_goa_folder(), "goa_parsed.txt", sep = "")
    write.table(goa, file = goa_annotations_dest, row.names = F, quote = F, sep = "\t")
    geneid2go <- topGO::readMappings(file = goa_annotations_dest)
    gene_names <- names(geneid2go)
    gene_list <- factor(as.integer(gene_names %in% gene_list))
    names(gene_list) <- gene_names
    godata <- new("topGOdata", ontology = ontology, allGenes = gene_list, annot = topGO::annFUN.gene2GO, gene2GO = geneid2go)

    # Performs enrichment test
    test_fisher <- new("classicCount", testStatistic = topGO::GOFisherTest, name = "Fisher test")
    result_fisher <- topGO::getSigGroups(godata, test_fisher)
    # Alternative: test.KS <- new('classicScore', testStatistic = GOKSTest, name = 'Kolmogorov-Smirnov test') Alternative: resultKS <- getSigGroups(godata, test.KS)

    # The subgraph induced by the top 5 GO terms identified by the elim algorithm for scoring GO terms for enrichment.  Rectangles indicate the 5 most significant terms. Rectangle color represents the
    # relative significance, ranging from dark red (most significant) to bright yellow (least significant). For each node, some basic information is displayed. The first two lines show the GO
    # identifier and a trimmed GO name. In the third line the raw p-value is shown. The forth line is showing the number of significant genes and the total number of genes annotated to the respective
    # GO term.  alertnative: showSigOfNodes(godata, score(result_fisher), firstSigNodes = 5, useInfo = 'all')

    results <- topGO::GenTable(godata, classicFisher = result_fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(names(result_fisher@score)))

    if (multiple_test_adjustment) {
        # Multiple test correction
        results$pvalue <- p.adjust(results$classicFisher, method = "fdr")
    } else {
        results$pvalue <- results$classicFisher
    }
    # removes not significant results (keeps those significant)
    results <- results[results$pvalue <= pvalue_threshold, ]

    futile.logger::flog.info("Finished computing GO enrichment")

    futile.logger::flog.info("Computes size of GO terms, based on Uniprot-GOA annotations")
    # Calculates frequency (size) of GO terms according to GOA (we use the whole uniprot to calculate size)
    sizes <- as.vector(sapply(as.character(results$GO.ID), FUN = get_goa_size))
    futile.logger::flog.info("Finished computing size of GO terms")

    # Returns results
    list(TopGOdata = godata, data.frame = data.frame(GO = as.character(results$GO.ID), pvalue = results$pvalue, name = results$Term, annotated_genes = results$Annotated, found_genes = results$Significant,
        expected_genes = results$Expected, size = sizes, row.names = NULL, stringsAsFactors = F))
}

# Calculates distances between GO terms based on the semantic similarity
get_go_distance_matrix <- function(enrichment_results, measure = "Wang", gene_list = NULL, ont = "BP", search_universe = "human") {
    # Creates a matrix with distances (i.e.: relations) between GO terms considering the genes associated between them
    go_distance_matrix <<- data.frame(matrix(0, nrow = length(enrichment_results$GO), ncol = length(enrichment_results$GO)), stringsAsFactors = F)
    rownames(go_distance_matrix) <- enrichment_results$GO
    colnames(go_distance_matrix) <- enrichment_results$GO
    if (measure %in% c("UI", "binary", "bray-curtis", "cosine")) {
        FUN <- get_functional_similarity
    } else {
        FUN <- get_semantic_similarity
    }

    # calculates similarity on pairwise comparisons
    pairwise_go_combinations <- combn(enrichment_results$GO, 2)
    apply(pairwise_go_combinations, 2, FUN = function(x) {
        go_distance_matrix[x[1], x[2]] <<- FUN(term1 = x[1], term2 = x[2], measure = measure, gene_list = gene_list, ont = ont, search_universe = search_universe)
    })

    # converts similarity matrix into distance matrix
    go_distance_matrix <- 1 - go_distance_matrix
    as.dist(t(go_distance_matrix))
    # we transpose the matrix as combn calculates similarities on the upper matrix
}

# Hierarchical clustering
hclust_clustering <- function(distance, distance_threshold) {
    # Clustering with hclust and a static distance threshold
    clustering <- hclust(distance)
    clusters <- cutree(clustering, h = distance_threshold)
    clusters
}

# PAM and silhouette clustering
pam_clustering <- function(distance) {
    # Clustering with PAM and silhouette
    max_clusters <- min(c(30, attr(distance, "Size") - 1))
    # checks if we have < 30 elements

    futile.logger::flog.info("Number of GO terms %d", attr(distance, "Size") - 1)
    futile.logger::flog.info("Looking best number of clusters between %d and %d", 2, max_clusters)
    sil_width <- numeric(max_clusters)
    for (k in 2:max_clusters) {
        sil_width[k] <- cluster::pam(distance, k)$silinfo$avg.width
    }
    k_best <- which.max(sil_width[2:max_clusters])
    clustering <- cluster::pam(distance, k_best)
    clusters <- clustering$clustering
    futile.logger::flog.info("Best number of clusters found of %d", k_best)
    clusters
}

# MDS
multidimensional_scaling <- function(distance) {
    number_elements <- attr(distance, "Size")
    futile.logger::flog.info("Applying MDS to %d elements", number_elements)
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

    list(x = x, y = y)
}

# Returns the semantic similarity score between the two terms
get_semantic_similarity <- function(term1, term2, measure = "Wang", gene_list = NULL, ont = "BP", search_universe = NULL) {
    res <- try(GOSemSim::mgoSim(term1, term2, organism = "human", measure = measure, ont = ont), silent = TRUE)
    if (class(res) == "try-error" | is.na(res)) {
        # Some terms are not found allegedly due to different versions of GO
        res <- 0
    }
    res
}

# Plots GO terms according to the enrichment results
cluster_and_plot <- function(enrichment_results, gene_list, topgo_data, method = "Wang", search_universe = "human", output_dir = ".", distance_threshold = 0.7, frequency_threshold = 0.05, clustering_method = "hclust",
    ont = "BP") {

    # Calculates the dissimilarity between GO terms
    if (method %in% c("Resnik", "Lin", "Rel", "Jiang", "Wang", "UI", "binary", "bray-curtis", "cosine", "all")) {

        # Prepares table
        enrichment_results <- cbind(enrichment_results, cluster = rep("", dim(enrichment_results)[1]), cluster_name = rep("", dim(enrichment_results)[1]))
        enrichment_results$cluster <- as.character(enrichment_results$cluster)
        enrichment_results$cluster_name <- as.character(enrichment_results$cluster_name)

        # Removes all terms with a frequency in GOA >5% enrichment_results = enrichment_results[enrichment_results$size <= frequency_threshold, ]

        if (method == "all") {
            methods <- c("Resnik", "Lin", "Rel", "Jiang", "Wang", "UI", "binary", "bray-curtis", "cosine")
        } else {
            methods <- c(method)
        }

        parent_output_dir <- output_dir
        for (methodn in methods) {

            output_dir <- paste(parent_output_dir, methodn, sep = "/")
            dir.create(output_dir)

            # Obtains the distance matrix between terms
            futile.logger::flog.info("Computing pairwise distance matrix on enriched GO terms with method %s", methodn)
            if (search_universe %in% c("human", "uniprot")) {
                go_distance <- get_go_distance_matrix(enrichment_results, measure = methodn, ont = ont, search_universe = search_universe)
            } else if (search_universe == "gene_list") {
                # similarity calculations are based on GOA associations limited to those genes in the gene list only applicable to binary, UI, Bray-Curtis and cosine distances
                go_distance <- get_go_distance_matrix(enrichment_results, measure = methodn, gene_list = gene_list, ont = ont, search_universe = search_universe)
            }
            futile.logger::flog.info("Finished computing distances")

            # Clustering
            futile.logger::flog.info("Computing clustering on enriched GO terms with method %s", clustering_method)
            if (clustering_method == "hclust") {
                clusters <- hclust_clustering(go_distance, distance_threshold)
            } else if (clustering_method == "pam") {
                clusters <- pam_clustering(go_distance)
            } else {
                stop("Clustering method not supported")
            }
            futile.logger::flog.info("Finished computing clustering")

            # Select cluster representatives
            futile.logger::flog.info("Selection of cluster representatives")
            for (cluster in unique(sort(clusters))) {
                group <- names(clusters)[clusters %in% cluster]
                # Selects as the cluster representative the most significantly associated term with a frequency <= 5% if any
                candidates <- enrichment_results[enrichment_results$GO %in% group & enrichment_results$size <= frequency_threshold, c("pvalue")]
                if (length(candidates) == 0) {
                  representative_idx <- which.min(enrichment_results[enrichment_results$GO %in% group, c("pvalue")])
                } else {
                  representative_idx <- which.min(enrichment_results[enrichment_results$GO %in% group & enrichment_results$size <= frequency_threshold, c("pvalue")])
                }
                enrichment_results[enrichment_results$GO %in% group, c("cluster")] <- group[representative_idx]
                enrichment_results[enrichment_results$GO %in% group, c("cluster_name")] <- enrichment_results[enrichment_results$GO == group[representative_idx], c("name")]
            }
            cluster_representatives <- enrichment_results[enrichment_results$GO == enrichment_results$cluster, ]

            # Multi Dimensional Scaling on all terms
            futile.logger::flog.info("Computing Multi Dimensional Scaling")
            mds_results <- multidimensional_scaling(go_distance)
            enrichment_results <- cbind(enrichment_results, mds_results$x)
            enrichment_results <- cbind(enrichment_results, mds_results$y)

            # Stores the results of clustering and MDS
            write.table(enrichment_results, file = paste(output_dir, "enrichment_results.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

            # Stores Euclidean distance to centroid to evaluate dispersion
            avg_x <- mean(mds_results$x)
            avg_y <- mean(mds_results$y)
            distances_to_centroid <- sqrt( ( mds_results$x - avg_x ) ^ 2 + ( mds_results$y - avg_y) ^ 2)
            write.table(data.frame(GO = names(distances_to_centroid), distance_to_centroid = distances_to_centroid), file = paste(output_dir, "distances_to_centroid.txt", sep = "/"), sep = "\t", quote = F,
                row.names = F)

            # Multi Dimensional Scaling just on clusters
            go_distance_m <- as.matrix(go_distance)
            cluster_distance <- as.dist(go_distance_m[rownames(go_distance_m) %in% cluster_representatives$GO, colnames(go_distance_m) %in% cluster_representatives$GO])
            mds_results <- multidimensional_scaling(cluster_distance)
            x <- mds_results$x
            y <- mds_results$y
            cluster_representatives <- cbind(cluster_representatives, x)
            cluster_representatives <- cbind(cluster_representatives, y)

            futile.logger::flog.info("Finished MDS")

            # Stores the results of clustering and MDS on representatives
            write.table(cluster_representatives, file = paste(output_dir, "cluster_representatives.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

            # Plots
            futile.logger::flog.info("Plots functional analysis results")
            try(plot_scatter(cluster_representatives, output_dir), silent = T)
            try(plot_table(cluster_representatives = cluster_representatives, output_dir), silent = T)
            try(plot_treemap(enrichment_results = enrichment_results, output_dir), silent = T)
            try(plot_graph(cluster_representatives = cluster_representatives, GOdata = topgo_data, output_dir), silent = T)
            futile.logger::flog.info("Finished plotting")
        }
    } else {
        # by default runs Wang semantic similarity
        stop("Non supported method.")
    }
}
