## Gets the folder to store GOA resources
get_goa_folder <- function() {
    get_package_folder("inst/GO/GOA/")
}

## Gets the folder to store GOA resources
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
        if (is.null(goa) || is.null(stored_ontology) || ontology != stored_ontology) {
            ## Downloads and stores as attribute human GOA
            futile.logger::flog.debug("Loading human GOA for the first time")
            uniKeys <- AnnotationDbi::keys(org.Hs.eg.db, keytype="SYMBOL")
            cols <- c("GOALL")
            goa_raw <- AnnotationDbi::select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="SYMBOL")
            goa_raw <- goa_raw[goa_raw$ONTOLOGYALL == ontology, ]
            goa_raw <- goa_raw[, c(1, 2)]
            colnames(goa_raw) <- c("Gene", "Term")
            goa <- new("GeneAnnotations", raw_annotations = goa_raw, name=paste("GOA-Human", ontology, sep=" "))
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
#' kegg = load_kegg()
load_kegg <- function() {

    kegg = NULL

    ## Loads stored attribute if any
    kegg <- attr(load_kegg, "human_kegg")
    if (is.null(kegg)) {
        ## Downloads and stores as attribute KEGG
        futile.logger::flog.debug("Loading human KEGG for the first time")
        uniKeys <- AnnotationDbi::keys(org.Hs.eg.db, keytype="SYMBOL")
        cols <- c("PATH")
        kegg_raw <- AnnotationDbi::select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="SYMBOL")
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
#' omim = load_omim()
load_omim <- function() {

    omim = NULL

    ## Loads stored attribute if any
    omim <- attr(load_omim, "omim")
    if (is.null(omim)) {
        ## Downloads and stores as attribute human OMIM
        futile.logger::flog.debug("Loading human OMIM for the first time")
        uniKeys <- AnnotationDbi::keys(org.Hs.eg.db, keytype="SYMBOL")
        cols <- c("OMIM")
        omim_raw <- AnnotationDbi::select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="SYMBOL")
        omim_raw <- omim_raw[, c(1, 2)]
        colnames(omim_raw) <- c("Gene", "Term")
        omim <- new("GeneAnnotations", raw_annotations = omim_raw, name="OMIM")
        attr(load_omim, "omim") <<- omim
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
        ## Downloads and stores as attribute HPO
        futile.logger::flog.debug("Loading HPO for the first time")
        file <- basename(hpo_annotations_url)
        file <- paste(hpo_annotations_folder, file, sep = "")
        download.file(hpo_annotations_url, file, quiet = TRUE)
        hpo_raw <- read.table(file, header = F, comment.char = "#", sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
        hpo_raw <- hpo_raw[, c(4, 1)]
        colnames(hpo_raw) <- c("Gene", "Term")
        ## Preprocess the HPO table, maps gene ids to entrez gene ids
        hpo <- new("GeneAnnotations", raw_annotations = hpo_raw, name="HPO")
        attr(load_hpo, "hpo") <<- hpo
    } else {
        futile.logger::flog.debug("Loading cached HPO")
    }

    return(hpo)
}

setClass("GeneAnnotations",
         representation(raw_annotations = "data.frame",
                        name = "character",
                        gene2term = "data.frame",
                        term2gene = "data.frame",
                        max_term_freq = "integer",
                        func_similarity_methods = "character"),
         prototype(
             raw_annotations = data.frame(),
             name = NA_character_),
         validity = function(object) {
             column_names = names(object@raw_annotations)
             stopifnot("Gene" %in% column_names & "Term" %in% column_names)
             }
         )

setMethod("initialize",
          signature(.Object = "GeneAnnotations"),
          function(.Object, raw_annotations, name){
              ## Removes entries with NA values
              .Object@raw_annotations <- raw_annotations[
                  !is.na(raw_annotations$Term) & !is.na(raw_annotations$Gene), ]
              .Object@name <- name
              .Object@gene2term <- aggregate(data = .Object@raw_annotations, Term ~ Gene, c)
              .Object@term2gene <- aggregate(data = .Object@raw_annotations, Gene ~ Term, c)
              .Object@max_term_freq <- max(sapply(.Object@term2gene$Gene, length))
              .Object@func_similarity_methods <- c("UI", "cosine", "bray-curtis", "binary")
              return(.Object)
          })


setGeneric("get_term_freq",
           signature = c("x", "term"),
           function(x, term) standardGeneric("get_term_freq"))

setMethod("get_term_freq", c("x" = "GeneAnnotations", "term" = "character"),
          function(x, term) {
              ## FIXME: this maximum size might need to be normalized as it is way too high
              ## somehow avoid considering terms very high in the hierarchy
              stopifnot(term %in% x@term2gene$Term)
              return(length(unlist(x@term2gene[x@term2gene$Term == term, ]$Gene)) / x@max_term_freq)
          })


setGeneric("get_functional_similarity",
           signature = c("x", "term1", "term2", "method"),
           function(x, term1, term2, method) standardGeneric("get_functional_similarity"))

setMethod("get_functional_similarity", c("x" = "GeneAnnotations",
                                         "term1" = "character",
                                         "term2" = "character",
                                         "method" = "character"),
          function(x, term1, term2, method) {
              #TODO: check valid methods
              stopifnot(term1 %in% x@term2gene$Term &&
                            term2 %in% x@term2gene$Term &&
                            method %in% x@func_similarity_methods)
              term1_genes <- x@term2gene[x@term2gene$Term == term1, c("Gene")][[1]]
              term2_genes <- x@term2gene[x@term2gene$Term == term2, c("Gene")][[1]]
              # calculates the union
              term_union <- union(term1_genes, term2_genes)
              if (length(term_union) == 0 | length(term1_genes) == 0 | length(term2_genes) == 0) {
                  # When no genes associated to any term they are disimilar
                  similarity <- 0
              } else {
                  # calculates the intersection
                  term_intersection <- intersect(term1_genes, term2_genes)
                  # calculates the different similarity metrics
                  if (method == "binary") {
                      term_xor <- term_union[!term_union %in% term_intersection]
                      similarity <- length(term_xor) / length(x@gene2term$Gene)
                  } else if (method == "UI") {
                      similarity <- length(term_intersection) / length(term_union)
                  } else if (method == "bray-curtis") {
                      similarity <- ( 2 * length(term_intersection) ) / ( length(term1_genes) + length(term2_genes) )
                  } else if (method == "cosine") {
                      similarity <- length(term_intersection) / ( sqrt(length(term1_genes) * length(term2_genes)) )
                  }
              }
              return(similarity)
          })
