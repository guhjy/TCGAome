library(biomaRt)

# Loads the Ensembl Biomart for Homo sapiens
load_biomart <- function ()
{
  if (! exists("ensembl"))
  {
    # Loads ensembl biomart for Homo sapiens
    #ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ensembl <<- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
    dataset="hsapiens_gene_ensembl"
    ensembl <<- useDataset(dataset, mart=ensembl)
  }
  
  ensembl
}

# Queries the biomart
query_biomart <- function (attributes, filters, values)
{
  # Example mapping HUGO gene names to Entrez gene IDs using BioMart
  #map_entrez_hgnc = getBM(attributes=c("entrezgene","hgnc_symbol"), filters=c("hgnc_symbol"), values=goa$Gene, mart=ensembl)
  getBM(attributes=attributes, filters=filters, values=values, mart=load_biomart())
}
