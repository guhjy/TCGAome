# TCGAome

## Running the tool
To get all available tumor types:
```
get_TCGAome_tumor_types()
```
Beware, that not all tumor types in the TCGA dataset are supported in TCGAome, only those having RNAseq and RPPA data. Also TCGAome works with the Firehose rdata release 20150402.

To run the tool:
```
run_TCGAome(tumor_types = c("BRCA", "OV"))
```

For non default configuration:

-param tumor_types Vector of tumor types to be analyzed from those available at TCGA (run RTCGAToolbox::getFirehoseDatasets() to see all available types)

-param run_pca_analysis Flag indicating if the PCA should run, this is part of the data pre-analysis (default: TRUE)

-param run_hclust_analysis Flag indicating if the hierarchical clustering analysis should run, this is part of the data pre-analysis (default: FALSE)

-param run_rgcca Flag indicating if the Regularized Generalized Canonical Correlation Analysis (RGCCA) should run (default: FALSE). Beware that this analysis is not intended for datasets with number of variables >> number of samples.

-param run_rcca Flag indicating if the Regularized Canonical Correlation Analysis (rCCA) should run (default: FALSE). Beware that this analysis is not intended for datasets with number of variables >> number of samples.

-param topN Indicates the top number of variables to select on MCIA and sPLS results (default: 5). It will select N variables on each of the data types, on each of the three first components and on each extreme of range, that is a maximum of 2*3*2*N, considering that there might be overlap between components.

-param spls_selection_method Indicates the method for variable selection on sPLS results. One of "correlation" or "loadings" (default: "loadings"). Loadings method will choose those variables maximizing variance across the samples, while correlation method will choose those variables with a higher correlation with other variables, that is those variables more distant to the origin in the correlation plot.

-param GO_similarity_measure The similarity measured employed to cluster and visualize GO terms, one in "Resnik", "Lin", "Rel", "Jiang", "Wang", "UI", "binary", "bray-curtis", "cosine", "all" (default: "Wang")

-param GOA_search_universe The GOA universe to use in functional similarity measures (i.e.: only aplicable to "binary", "bray-curtis", "cosine", "UI"), one in "human", "uniprot", "gene_list" (default: "human")

-param enrichment_significance_threshold The threshold to consider a GO enrichment result as significant (default: 0.01)

-param GO_ontology The GO ontology on which to perform the enrichment, one of "BP", "MF", "CC" (default: BP)

-param multiple_test_adjustment Flag to run multiple test adjustment on GO enrichment results (default: F)


## Installation

On Linux systems you might need additional installation of some Cairo and RGL libs for the graphical output. On Ubuntu:
```
apt-get install libcairo2-dev
apt-get install libxt-dev
apt-get install libglu1-mesa-dev
```

Many of TCGAome dependencies are not in the CRAN repository and they might require a previous manual installation from Bioconductor:
```
source("http://bioconductor.org/biocLite.R")
biocLite("omicade4")
biocLite("mixOmics")
biocLite("RGCCA")
biocLite("GOSemSim")
biocLite("RTCGAToolbox")
biocLite("ReactomePA")
biocLite("biomaRt")
biocLite("topGO")
```

Installs the R development tools:
```
install.packages(devtools)
library(devtools)
```

Installs the ggbiplot package from GitHub:
```
install_github("ggbiplot", "vqv")
```

Install the TCGAome package from GitHub:
```
install_github('TCGAome', 'priesgo')
```

