

setClass("GeneAnnotations",
         representation(raw_annotations = "data.frame",
                        name = "character",
                        gene2term = "data.frame",
                        term2gene = "data.frame",
                        max_term_ic = "integer"),
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
              .Object@raw_annotations <- raw_annotations
              .Object@name <- name
              .Object@gene2term <- aggregate(data = .Object@raw_annotations, Term ~ Gene, c)
              .Object@term2gene <- aggregate(data = .Object@raw_annotations, Gene ~ Term, c)
              .Object@max_term_ic <- max(sapply(.Object@term2gene$Gene, length))
              return(.Object)
          })


setGeneric("get_term_size",
           signature = c("x", "term"),
           function(x, term) standardGeneric("get_term_size"))

setMethod("get_term_size", c("x" = "GeneAnnotations", "term" = "character"),
          function(x, term) {
              ## FIXME: this maximum size might need to be normalized as it is way too high
              ## somehow avoid considering terms very high in the hierarchy
              return(length(unique(x@term2gene[x@term2gene$Term == term, ]$Gene[[1]])) / x@max_term_ic)
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
                      similarity <- length(term_xor) / length(unique(goa$Gene))
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
