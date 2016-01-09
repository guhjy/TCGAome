library(futile.logger)

#!/usr/bin/env Rscript
main <- function(){
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
  # to load from RStudio
  # script.basename = getwd()

  # Loads all code in this folder
  #for (nm in list.files(script.basename, pattern = "\\.[RrSsQq]$")) {
  for (nm in list.files("./R/", pattern = "\\.[RrSsQq]$")) {
    if (nm != "main.r"){
      source(file.path("./R", nm))
    }
  }


  # Sets logging
  flog.threshold(DEBUG)
  flog.appender(appender.file("../logs/log"), name="log")

  # Runs the pipeline
  #pipeline(c("BRCA", "OV"))
}


configure_logging <- function(folder){
  flog.threshold(DEBUG)
  flog.appender(appender.file(paste(folder,"TCGAome.log", sep="/")), name="log")
}


get_package_folder <- function(relative_path){
  paste(system.file(package = "TCGAome"), relative_path, sep="/")
}

get_results_folder <- function()
{
  timestamp = format(Sys.time(),"%Y%m%d_%H%M%S")
  results_folder = get_package_folder(paste("results", timestamp, sep="_"))
  dir.create(results_folder)

  results_folder
}

# Bioconductor packages are not loaded by import statements, so this functions loads them explicitly.
loads_bioc_packages <- function()
{
  library(ReactomePA)
  library(RTCGAToolbox)
  library(mixOmics)
  library(omicade4)
  library(RGCCA)
  library(biomaRt)
  library(ontoCAT)
  library(topGO)
  library(GOSemSim)
}


