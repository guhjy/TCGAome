## Gets the folder to store GOA resources
get_goa_folder <- function() {
    get_package_folder("inst/GO/GOA/")
}

## Gets the folder to store GOA resources
get_hpo_folder <- function() {
    get_package_folder("inst/HPO/")
}


#' Reads Gene Ontology Annotations (GOA) for Human genes and for all Uniprot.
#' No electronically inferred annotations (IEA)
#' @param search_universe Two possible values "human" or "uniprot".
#' The first provides associations to GO of human genes,
#' while the second considers all genes registered in Uniprot.
#' [default: "human"]
#' @param ontology The ontology within Gene Ontology: "BP", "CC" or "MF".
#' [default: "BP"]
#' @param goa_uniprot_annotations_url The URL to Uniprot GOA.
#' [default: "http://geneontology.org/gene-associations/gene_association.goa_uniprot_noiea.gz"]
#'
#' @keywords TCGAome
#' @export
#' @examples
#' goa = load_goa()
load_goa <- function(
    search_universe = "human",
    ontology = "BP",
    goa_uniprot_annotations_url =
        "http://geneontology.org/gene-associations/gene_association.goa_uniprot_noiea.gz") {

    if (! search_universe %in% c("human", "uniprot")) {
        stop("Invalid value for search_universe. Expected one of ['human', 'uniprot']")
    }
    if (! ontology %in% c("BP", "CC", "MF")) {
        stop("Invalid value for ontology. Expected one of ['BP', 'CC', 'MF']")
    }

    ## get and create working folder
    goa_annotations_folder <- TCGAome::get_goa_folder()
    goa = NULL

    if (search_universe == "human") {
        ## Downloads and stores as attribute human GOA
        futile.logger::flog.info("Loading human GOA...")
        uniKeys <- AnnotationDbi::keys(
            org.Hs.eg.db::org.Hs.eg.db,
            keytype="SYMBOL")
        cols <- c("GOALL")
        goa_raw <- AnnotationDbi::select(
            org.Hs.eg.db::org.Hs.eg.db,
            keys=uniKeys,
            columns=cols,
            keytype="SYMBOL")
        goa_raw <- goa_raw[goa_raw$ONTOLOGYALL == ontology, ]
        goa_raw <- goa_raw[, c(1, 2)]
        colnames(goa_raw) <- c("Gene", "Term")
        name <- paste("GOA-Human", ontology, sep=" ")
    } else if (search_universe == "uniprot") {
        ## Downloads and stores as attribute uniprot GOA
        futile.logger::flog.info("Loading uniprot GOA...")
        file <- basename(goa_uniprot_annotations_url)
        file <- paste(goa_annotations_folder, file, sep = "")
        download.file(goa_uniprot_annotations_url, file, quiet = TRUE)
        goa_raw <- read.table(
            file, header = F, comment.char = "!",
            sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
        goa_raw <- goa_raw[, c(3, 5)]
        colnames(goa_raw) <- c("Gene", "Term")
        name <- "GOA-Uniprot"
    }
    goa <- new("GeneAnnotations",
               raw_annotations = goa_raw,
               name = name)
    futile.logger::flog.info(
        "Loaded %d GO terms and %d genes",
        nrow(goa@term2gene),
        nrow(goa@gene2term))
    return(goa)
}

#' Reads KEGG pathways for Human genes and for all Uniprot.
#' No electronically inferred annotations (IEA)
#'
#' @keywords TCGAome
#' @export
#' @examples
#' kegg = load_kegg()
load_kegg <- function() {
    ## Downloads and stores as attribute KEGG
    futile.logger::flog.info("Loading human KEGG...")
    uniKeys <- AnnotationDbi::keys(
        org.Hs.eg.db::org.Hs.eg.db,
        keytype="SYMBOL")
    cols <- c("PATH")
    kegg_raw <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys=uniKeys,
        columns=cols,
        keytype="SYMBOL")
    kegg_raw <- kegg_raw[, c(1, 2)]
    colnames(kegg_raw) <- c("Gene", "Term")
    kegg <- new("GeneAnnotations",
                raw_annotations = kegg_raw,
                name="KEGG-Human")
    futile.logger::flog.info(
        "Loaded %d KEGG pathways and %d genes",
        nrow(kegg@term2gene),
        nrow(kegg@gene2term))
    return(kegg)
}

#' Reads OMIM diseases for Human genes.
#'
#' @keywords TCGAome
#' @export
#' @examples
#' omim = load_omim()
load_omim <- function() {
    ## Downloads and stores as attribute human OMIM
    futile.logger::flog.info("Loading human OMIM...")
    uniKeys <- AnnotationDbi::keys(
        org.Hs.eg.db::org.Hs.eg.db,
        keytype="SYMBOL")
    cols <- c("OMIM")
    omim_raw <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys=uniKeys,
        columns=cols,
        keytype="SYMBOL")
    omim_raw <- omim_raw[, c(1, 2)]
    colnames(omim_raw) <- c("Gene", "Term")
    omim <- new("GeneAnnotations",
                raw_annotations = omim_raw,
                name="OMIM")
    futile.logger::flog.info(
        "Loaded %d OMIM diseases and %d genes",
        nrow(omim@term2gene),
        nrow(omim@gene2term))
    return(omim)
}

#' Reads Human Phenotype Ontology (HPO) annotations.
#' @param hpo_annotations_url The URL to HPO annotations to genes.
#' [default: "http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_TYPICAL_FEATURES_phenotype_to_genes.txt"]
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

    ## Downloads and stores as attribute HPO
    futile.logger::flog.info("Loading HPO...")
    file <- basename(hpo_annotations_url)
    file <- paste(hpo_annotations_folder, file, sep = "")
    download.file(hpo_annotations_url, file, quiet = TRUE)
    hpo_raw <- read.table(
        file, header = F, comment.char = "#",
        sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
    hpo_raw <- hpo_raw[, c(4, 1)]
    colnames(hpo_raw) <- c("Gene", "Term")
    ## Preprocess the HPO table, maps gene ids to entrez gene ids
    hpo <- new("GeneAnnotations",
               raw_annotations = hpo_raw,
               name="HPO")
    futile.logger::flog.info(
        "Loaded %d HPO phenotypes and %d genes",
        nrow(hpo@term2gene),
        nrow(hpo@gene2term))
    return(hpo)
}

#' GeneAnnotations class
#'
#' An S4 class to represent a set of gene annotations built from
#' tall-skinny table with the columns "Gene" and "Term". This class provides
#' the functionality to compute several useful metrics based on the gene
#' annotations: distance matrix, term frequency of annotation,
#' functional similarity and enrichment.
#'
#' @slot raw_annotations A data.frame having columns "Gene" and "Term" that
#' contains all annotations.
#' @slot name The annotations name
#' @slot gene2term A data.frame that stores the pivoted annotations across genes
#' @slot term2gene A data.frame that stores the pivoted annotations across terms
#' @slot max_term_annotations The maximum frequency of annotations, that is the maximum
#' number of terms associated to any gene
#' @slot func_similarity_methods Stores the supported similarity methods
#'
#' @export
#'
setClass("GeneAnnotations",
         representation(raw_annotations = "data.frame",
                        name = "character",
                        gene2term = "data.frame",
                        term2gene = "data.frame",
                        max_term_annotations = "integer",
                        func_similarity_methods = "character"),
         #prototype(
         #    raw_annotations = data.frame(),
         #    name = NA_character_),
         validity = function(object) {
                 errors <- character()
                 # Check that the annotations data frame has the expected column names
                 column_names = names(object@raw_annotations)
                 if (! "Gene" %in% column_names | ! "Term" %in% column_names){
                     msg <- paste("Column names should contain columns 'Gene' and 'Term'. Columns found [", str(column_names), "]", sep = "")
                     errors <- c(errors, msg)
                 }
                 # Check that the annotations data frame is not empty
                 if (dim(object@raw_annotations)[1] == 0){
                     msg <- paste("Empty annotations")
                     errors <- c(errors, msg)
                 }
                 # Check that there are at least 2 terms
                 if (length(unique(object@raw_annotations$Term)) < 2){
                     msg <- paste("Annotations with at least 2 terms are required")
                     errors <- c(errors, msg)
                 }

                 if (length(errors) == 0) TRUE else errors
             }
         )

#' GeneAnnotations constructor
#'
#' The constructor for GeneAnnotations
#'
#' @param raw_annotations A data.frame having columns "Gene" and "Term" that
#' contains all annotations.
#' @param name The annotations name
#'
#' @return The GeneAnnotations object
#'
#' @export
#'
#' @examples
#' uniKeys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype="SYMBOL")
#' cols <- c("PATH")
#' kegg_raw <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="SYMBOL")
#' kegg_raw <- kegg_raw[, c(1, 2)]
#' colnames(kegg_raw) <- c("Gene", "Term")
#' kegg <- GeneAnnotations(raw_annotations = kegg_raw, name="KEGG-Human")
GeneAnnotations <- function(...) new("GeneAnnotations",...)

setMethod("initialize",
          signature(.Object = "GeneAnnotations"),
          function(.Object, raw_annotations, name){
              ## Avoids factor columns
              raw_annotations$Gene <- as.character(raw_annotations$Gene)
              raw_annotations$Term <- as.character(raw_annotations$Term)
              ## Removes entries with NA values
              .Object@raw_annotations <- raw_annotations[
                  !is.na(raw_annotations[, "Term"]) & !is.na(raw_annotations[, "Gene"]), ]
              ## Sets slots
              .Object@name <- name
              .Object@gene2term <- aggregate(data = .Object@raw_annotations,
                            Term ~ Gene, c)
              rownames(.Object@gene2term) <- .Object@gene2term$Gene
              .Object@term2gene <- aggregate(data = .Object@raw_annotations,
                            Gene ~ Term, c)
              rownames(.Object@term2gene) <- .Object@term2gene[, c("Term")]
              .Object@max_term_annotations <- max(sapply(.Object@term2gene[, c("Gene")], length))
              .Object@func_similarity_methods <- c("UI", "cosine", "bray-curtis", "binary")

              ## Checks object validity
              validObject(.Object)

              return(.Object)
          })


#' get_term_freq()
#'
#' Function that gets the frequency of annotation of a specific term. That is
#' the number of genes this term is associated to divided by the maximum term
#' annotation.
#'
#' @param x The GeneAnnotations class on which the method will run.
#' @param term The term on which we want the annotation frequency.

#' @return The term's annotations frequency
#'
#' @export
#'
#' @examples
#' kegg <- TCGAome::load_kegg()
#' random_kegg_term = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' get_term_freq(kegg, random_kegg_term)
setGeneric("get_term_freq",
           signature = c("x", "term"),
           function(x, term) standardGeneric("get_term_freq"))

#' @aliases get_term_freq
#' @export
setMethod("get_term_freq", c("x" = "GeneAnnotations", "term" = "character"),
          function(x, term) {
              ## FIXME: this maximum size might need to be normalized as it is way too high
              ## somehow avoid considering terms very high in the hierarchy
              #if (! term %in% x@term2gene$Term) {
              #    stop(paste("Term ", term, " not present in the annotations", sep = ""))
              #}
              return(length(unlist(
                  x@term2gene[term, ]$Gene))
                  / x@max_term_annotations)
          })

#' get_functional_similarity()
#'
#' Function that gets the functional similarity between two terms. That is
#' a measure of the shared terms in the annotations resources. This functional
#' similarity is implemented as a distance between binary vectors.
#'
#' @param x The GeneAnnotations class on which the method will run.
#' @param term1 The first term on which to calculate the functional similarity.
#' @param term2 The second term on which to calculate the functional similarity.
#' @param distance_measure The binary distance method (one of UI, binary,
#' bray-curtis, cosine)

#' @return The functional similarity between term1 and term2
#'
#' @export
#'
#' @examples
#' kegg <- TCGAome::load_kegg()
#' random_kegg_term1 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' random_kegg_term2 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' get_functional_similarity(kegg, random_kegg_term1, random_kegg_term2, "cosine")
#' get_functional_similarity(kegg, random_kegg_term1, random_kegg_term2, "binary")
#' get_functional_similarity(kegg, random_kegg_term1, random_kegg_term2, "UI")
#' get_functional_similarity(kegg, random_kegg_term1, random_kegg_term2, "bray-curtis")
setGeneric("get_functional_similarity",
           signature = c("x", "term1", "term2", "distance_measure"),
           function(x, term1, term2, distance_measure) standardGeneric("get_functional_similarity"))

#' @aliases get_functional_similarity
#' @export
setMethod("get_functional_similarity", c("x" = "GeneAnnotations",
                                         "term1" = "character",
                                         "term2" = "character",
                                         "distance_measure" = "character"),
          function(x, term1, term2, distance_measure) {
              if (distance_measure == "binary") {
                  similarity <- TCGAome::get_binary_similarity(x, term1, term2)
              } else if (distance_measure == "UI") {
                  similarity <- TCGAome::get_ui_similarity(x, term1, term2)
              } else if (distance_measure == "bray-curtis") {
                  similarity <- TCGAome::get_bc_similarity(x, term1, term2)
              } else if (distance_measure == "cosine") {
                  similarity <- TCGAome::get_cosine_similarity(x, term1, term2)
              } else {
                  stop(paste("Non supported distance measure ", distance_measure, sep=""))
              }
              return(similarity)
          })


#' get_binary_similarity()
#'
#' Function that gets the binary functional similarity between two terms. That is
#' a measure of the shared terms in the annotations resources.
#'
#' @param x The GeneAnnotations class on which the method will run.
#' @param term1 The first term on which to calculate the functional similarity.
#' @param term2 The second term on which to calculate the functional similarity.

#' @return The binary functional similarity between term1 and term2
#'
#' @export
#'
#' @examples
#' kegg <- TCGAome::load_kegg()
#' random_kegg_term1 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' random_kegg_term2 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' get_binary_similarity(kegg, random_kegg_term1, random_kegg_term2)
#setGeneric("get_binary_similarity",
#           signature = c("x", "term1", "term2"),
#           function(x, term1, term2) standardGeneric("get_binary_similarity"))

#' @aliases get_binary_similarity
#' @export
#setMethod("get_binary_similarity", c("x" = "GeneAnnotations",
#                                         "term1" = "character",
#                                         "term2" = "character"),
get_binary_similarity <- function(x, term1, term2) {
    column_names <- names(x@term2gene)
    term1_genes <- x@term2gene[term1, ]$Gene[[1]]
    term2_genes <- x@term2gene[term2, ]$Gene[[1]]
    if (length(term1_genes) == 0 | length(term2_genes) == 0) {
        # When no genes associated to any term they are disimilar
        similarity <- 0.0
    } else {
        term_xor <- c(term1_genes[is.na(fmatch(term1_genes, term2_genes))],
                      term2_genes[is.na(fmatch(term2_genes, term1_genes))])
        similarity <- length(term_xor) / length(x@gene2term$Gene)
    }
    return(similarity)
}

get_binary_similarity_bc <- compiler::cmpfun(get_binary_similarity)


#' get_ui_similarity()
#'
#' Function that gets the Union-Intersection functional similarity between two terms. That is
#' a measure of the shared terms in the annotations resources.
#'
#' @param x The GeneAnnotations class on which the method will run.
#' @param term1 The first term on which to calculate the functional similarity.
#' @param term2 The second term on which to calculate the functional similarity.

#' @return The UI functional similarity between term1 and term2
#'
#' @export
#'
#' @examples
#' kegg <- TCGAome::load_kegg()
#' random_kegg_term1 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' random_kegg_term2 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' get_ui_similarity(kegg, random_kegg_term1, random_kegg_term2)
#setGeneric("get_ui_similarity",
#           signature = c("x", "term1", "term2"),
#           function(x, term1, term2) standardGeneric("get_ui_similarity"))

#' @aliases get_ui_similarity
#' @export
#setMethod("get_ui_similarity", c("x" = "GeneAnnotations",
#                                     "term1" = "character",
#                                     "term2" = "character"),
get_ui_similarity <- function(x, term1, term2) {
    column_names <- names(x@term2gene)
    term1_genes <- x@term2gene[term1, ]$Gene[[1]]
    term2_genes <- x@term2gene[term2, ]$Gene[[1]]
    if (length(term1_genes) == 0 | length(term2_genes) == 0) {
        # When no genes associated to any term they are disimilar
        similarity <- 0.0
    } else {
        term_intersection <- term1_genes[!is.na(fmatch(term1_genes, term2_genes))]
        term_union <- c(term1_genes,
                        term2_genes[is.na(fmatch(term2_genes, term1_genes))])
        similarity <- length(term_intersection) / length(term_union)
    }
    return(similarity)
}

get_ui_similarity_bc <- compiler::cmpfun(get_ui_similarity)

#' get_bc_similarity()
#'
#' Function that gets the Bray-Curtis functional similarity between two terms. That is
#' a measure of the shared terms in the annotations resources.
#'
#' @param x The GeneAnnotations class on which the method will run.
#' @param term1 The first term on which to calculate the functional similarity.
#' @param term2 The second term on which to calculate the functional similarity.

#' @return The Bray-Curtis functional similarity between term1 and term2
#'
#' @export
#'
#' @examples
#' kegg <- TCGAome::load_kegg()
#' random_kegg_term1 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' random_kegg_term2 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' get_bc_similarity(kegg, random_kegg_term1, random_kegg_term2)
#setGeneric("get_bc_similarity",
#           signature = c("x", "term1", "term2"),
#           function(x, term1, term2) standardGeneric("get_bc_similarity"))

#' @aliases get_bc_similarity
#' @export
#setMethod("get_bc_similarity", c("x" = "GeneAnnotations",
#                                 "term1" = "character",
#                                 "term2" = "character"),
get_bc_similarity <- function(x, term1, term2) {
    column_names <- names(x@term2gene)
    term1_genes <- x@term2gene[term1, ]$Gene[[1]]
    term2_genes <- x@term2gene[term2, ]$Gene[[1]]
    if (length(term1_genes) == 0 | length(term2_genes) == 0) {
        # When no genes associated to any term they are disimilar
        similarity <- 0.0
    } else {
        term_intersection <- term1_genes[!is.na(fmatch(term1_genes, term2_genes))]
        similarity <- ( 2 * length(term_intersection) ) / ( length(term1_genes) + length(term2_genes) )
    }
    return(similarity)
}

get_bc_similarity_bc <- compiler::cmpfun(get_bc_similarity)

#' get_cosine_similarity()
#'
#' Function that gets the cosine functional similarity between two terms. That is
#' a measure of the shared terms in the annotations resources.
#'
#' @param x The GeneAnnotations class on which the method will run.
#' @param term1 The first term on which to calculate the functional similarity.
#' @param term2 The second term on which to calculate the functional similarity.

#' @return The cosine functional similarity between term1 and term2
#'
#' @export
#'
#' @examples
#' kegg <- TCGAome::load_kegg()
#' random_kegg_term1 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' random_kegg_term2 = kegg@term2gene$Term[runif(1, max = length(kegg@term2gene$Term))]
#' get_cosine_similarity(kegg, random_kegg_term1, random_kegg_term2)
get_cosine_similarity <- function(x, term1, term2) {
    column_names <- names(x@term2gene)
    term1_genes <- x@term2gene[term1, ]$Gene[[1]]
    term2_genes <- x@term2gene[term2, ]$Gene[[1]]
    if (length(term1_genes) == 0 | length(term2_genes) == 0) {
        # When no genes associated to any term they are disimilar
        similarity <- 0.0
    } else {
        term_intersection <- term1_genes[!is.na(fmatch(term1_genes, term2_genes))]
        similarity <- length(term_intersection) /  sqrt(length(term1_genes) * length(term2_genes))
    }
    return(similarity)
}

get_cosine_similarity_bc <- compiler::cmpfun(get_cosine_similarity)


.createCluster = function(noCores, logfile = "/dev/null", export = NULL, lib = NULL) {
    cl <- parallel::makeCluster(noCores, type = "SOCK", outfile = logfile)
    if(!is.null(export)) clusterExport(cl, export)
    if(!is.null(lib)) {
        l_ply(lib, function(dum) {
            parallel::clusterExport(cl, "dum", envir = environment())
            parallel::clusterEvalQ(cl, library(dum, character.only = TRUE))
        })
    }
    doSNOW::registerDoSNOW(cl)
    return(cl)
}

#' get_term_distance_matrix()
#'
#' Function that returns the distance matrix between all the terms.
#'
#' @param x The GeneAnnotations class on which the method will run.
#' @param distance_measure The binary distance method (one of UI, binary,
#' bray-curtis, cosine)
#' @param terms_subset A subset of terms on which to calculate the distance matrix

#' @return The distance matrix for the subset if provided or for all the terms
#' otherwise
#'
#' @export
#'
#' @examples
#' kegg <- TCGAome::load_kegg()
#' get_term_distance_matrix(kegg, "cosine")
setGeneric("get_term_distance_matrix",
           signature = c("x", "distance_measure"),
           function(x, distance_measure, ..., terms_subset = NULL) standardGeneric("get_term_distance_matrix"))

#' @aliases get_term_distance_matrix
#' @export
setMethod("get_term_distance_matrix", c("x" = "GeneAnnotations",
                                        "distance_measure" = "character"),
          function(x, distance_measure, ..., terms_subset = NULL) {

              if (! distance_measure %in% x@func_similarity_methods) {
                  stop(paste("Non supported distance measure ", distance_measure, sep=""))
              }

              if (!exists("terms_subset") || is.null(terms_subset)) {
                  warning("Compute distance matrix on all terms. This may be time consuming")
                  terms_subset <- x@term2gene$Term
              }

              if (length(terms_subset) < 2) {
                  stop("Cannot compute distance matrix on less than two terms")
              }

              ## Initializes to 1 as diagonal elements will not be evaluated.
              similarity_matrix <- data.frame(matrix(1,
                                                     nrow = length(terms_subset),
                                                     ncol = length(terms_subset)),
                                              stringsAsFactors = F)
              rownames(similarity_matrix) <- terms_subset
              colnames(similarity_matrix) <- terms_subset

              # Calculates similarity on pairwise comparisons without repetition
              pairwise_term_combn <- combn(terms_subset, 2)

              # Set the functional similarity function
              if (distance_measure == "binary") {
                  get_similarity <- get_binary_similarity_bc
              } else if (distance_measure == "UI") {
                  get_similarity <- get_ui_similarity_bc
              } else if (distance_measure == "bray-curtis") {
                  get_similarity <- get_bc_similarity_bc
              } else if (distance_measure == "cosine") {
                  get_similarity <- get_cosine_similarity_bc
              }
              require(fastmatch)

              # Computes the distance
              #nodes <- parallel::detectCores()
              #cl <- .createCluster(nodes)
              #tallskinny_dist <- plyr::aaply(.data = pairwise_term_combn, .margin = 2, .fun = function(y) {
              #    get_functional_similarity(
              #        x, term1 = y[1], term2 = y[2], distance_measure = distance_measure)
              #    }, .parallel = TRUE)
              #stopCluster(cl)
              #
              # Snow
              #nodes <- parallel::detectCores()
              #cl <- snow::makeCluster(4, type = "SOCK")
              #tallskinny_dist <- snow::parApply(
              #    cl, pairwise_term_combn, 2,
              #    FUN = function(y) {
              #        get_similarity(x, term1 = y[1], term2 = y[2])
              #    })
              #snow::stopCluster(cl)
              #
              tallskinny_dist <- apply(
                  pairwise_term_combn, 2,
                  FUN = function(y) {
                      get_similarity(x, term1 = y[1], term2 = y[2])
                      })
              idx <- rbind(cbind(fmatch(pairwise_term_combn[1 , ], terms_subset),
                                 fmatch(pairwise_term_combn[2 , ], terms_subset)),
                           cbind(fmatch(pairwise_term_combn[2 , ], terms_subset),
                                 fmatch(pairwise_term_combn[1 , ], terms_subset)))
              similarity_matrix[idx] <- tallskinny_dist

              # Converts similarity matrix into distance matrix
              distance_matrix <- 1 - similarity_matrix
              return(as.dist(distance_matrix))
          })

#' get_enrichment()
#'
#' Function that returns an object of class GeneListEnrichment with the
#' enrichment results on a given gene list on this gene annotations.
#'
#' @param x The GeneAnnotations class on which the method will run.
#' @param gene_list The list of genes on which to calculate the enrichment.

#' @return The enrichment results stored in an object of class
#' GeneListEnrichment
#'
#' @export
#'
#' @examples
#' kegg <- TCGAome::load_kegg()
#' gene_list <- c("ZNF638", "HNRNPU", "PPIAL4G", "RAPH1", "USP7", "SUMO1P3",
#' "TMEM189.UBE2V1", "ZNF837", "LPCAT4", "ZFPL1", "STAT3", "XRCC1",
#' "STMN1", "PGR", "RB1", "KDR", "YBX1", "YAP1", "FOXO3", "SYK", "RAB17",
#' "TTC8", "SLC22A5", "C3orf18", "ANKRA2", "LBR", "B3GNT5", "ANP32E",
#' "JOSD1", "ZNF695", "ESR1", "INPP4B", "PDK1", "TSC2", "AR", "HSPA1A",
#' "CDH3", "SMAD4", "CASP7", "GMPS", "NDC80", "EZH2", "MELK", "CDC45",
#' "CRY2", "KLHDC1", "MEIS3P1", "FBXL5", "EHD2", "CCNB1", "GSK3A",
#' "DVL3", "NFKB1", "COL6A1", "CCND1", "BAK1")
#' get_enrichment(kegg, gene_list)
setGeneric("get_enrichment",
           signature = c("x", "gene_list"),
           function(x, gene_list) standardGeneric("get_enrichment"))

#' @aliases get_enrichment
#' @export
setMethod("get_enrichment", c("x" = "GeneAnnotations", "gene_list" = "character"),
          function(x, gene_list) {
              return(TCGAome::GeneListEnrichment(gene_list = gene_list, gene_annotations = x))
          })

#' plot()
#'
#' Plots the gene annotations
#'
#' @param x The GeneAnnotations class on which the method will run.

#' @return The ggplot2 plot
#'
#' @export
#'
#' @examples
#' TCGAome::plot(hpo)
setGeneric("plot",
           signature = c("x"),
           function(x) standardGeneric("plot"))

#' @aliases plot
#' @export
setMethod("plot",
          c("x" = "GeneAnnotations"),
          function(x) {

              gene_dist <- as.data.frame(
                  table(
                      unlist(lapply(x@gene2term$Term, FUN = length))),
                  stringsAsFactors = F)
              gene_dist$class <- "gene"
              names(gene_dist) <- c("nterms", "ngenes", "class")

              term_dist <- as.data.frame(
                  table(
                      unlist(lapply(x@term2gene$Gene, FUN = length))),
                  stringsAsFactors = F)
              term_dist$class <- "term"
              names(term_dist) <- c("ngenes", "nterms", "class")

              annotation_dist <- rbind(gene_dist, term_dist)
              annotation_dist$ngenes <- as.numeric(annotation_dist$ngenes)
              annotation_dist$nterms <- as.numeric(annotation_dist$nterms)

              gene_dist <- data.frame(
                  class = "gene",
                  associations = unlist(lapply(x@gene2term$Term, FUN = length))
              )
              term_dist <- data.frame(
                  class = "term",
                  associations = unlist(lapply(x@term2gene$Gene, FUN = length))
              )
              annotation_dist2 <- rbind(gene_dist, term_dist)

              ## Plots gene associations
              mb <- as.numeric(1:10 %o% 10 ^ (0:4))
              plot1 <- ggplot2::ggplot(data = annotation_dist,
                                      ggplot2::aes(x = nterms, y = ngenes)) +
                  ggplot2::geom_point(ggplot2::aes(shape = class,
                                                   color = class)) +
                  ggplot2::scale_y_log10(
                      breaks = c(1, 10, 100, 1000, 10000),
                      minor_breaks = mb) +
                  ggplot2::scale_x_log10(
                      breaks = c(1, 10, 100, 1000, 10000),
                      minor_breaks = mb) +
                  ggplot2::scale_shape_manual(values=c(1, 2)) +
                  ggplot2::labs(x = "# terms",
                                y = "# genes",
                                colour = NULL, shape = NULL,
                                title = paste(x@name, "joint dist", sep = " - ")) +
                  ggplot2::scale_color_brewer(
                      palette = "Set1") +
                  ggplot2::theme_bw()

              plot2 <- ggplot2::ggplot(data = annotation_dist2,
                                      ggplot2::aes(class, associations, colour = class)) +
                  ggplot2::geom_violin() +
                  ggplot2::scale_y_log10(
                      breaks = c(1, 10, 100, 1000, 10000),
                      minor_breaks = mb) +
                  ggplot2::labs(x = NULL,
                                y = "Associations",
                                colour = NULL, shape = NULL,
                                title = paste(x@name, "disjoint dist", sep = " - ")) +
                  ggplot2::scale_color_brewer(
                      palette = "Set1") +
                  ggplot2::theme_bw()

              cowplot::plot_grid(plot1, plot2, nrow = 2)
          })
