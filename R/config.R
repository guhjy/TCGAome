

# Configures loggers to be used
configure_logging <- function(folder){
  flog.threshold(DEBUG)
  flog.appender(appender.tee(paste(folder,"TCGAome.log", sep="/")))
}


get_package_folder <- function(relative_path){
  paste(gsub("\\\\", "/", base::system.file(package = "TCGAome")), relative_path, sep="/")
}

get_results_folder <- function()
{
  timestamp = format(Sys.time(),"%Y%m%d_%H%M%S")
  results_folder = get_package_folder(paste("results", timestamp, sep="_"))
  dir.create(results_folder)

  results_folder
}

# Bioconductor packages are not loaded by import statements, so this functions loads them explicitly.
loads_dependencies <- function()
{
  library(ReactomePA)
  library(RTCGAToolbox)
  library(mixOmics)
  library(omicade4)
  library(biomaRt)
  library(ontoCAT)
  library(topGO)
  library(GOSemSim)
  library(treemap)
  library(gridExtra)
  library(ggplot2)
  library(scales)
  library(grid)
  library(VennDiagram)
  library(ggbiplot)
  library(gplots)
  library(Cairo)
  library(reshape)
  library(cluster)
  library(dplyr)
  library(futile.logger)
}


