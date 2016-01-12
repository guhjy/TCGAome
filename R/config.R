

# Configures loggers to be used
configure_logging <- function(folder){
  flog.threshold(DEBUG)
  flog.appender(appender.tee(paste(folder,"TCGAome.log", sep="/")))
}

# Gets the absolute path for the package folder
get_package_folder <- function(relative_path){
  paste(gsub("\\\\", "/", base::system.file(package = "TCGAome")), relative_path, sep="/")
}

# Gets the results folder with a unique timestamp
get_results_folder <- function()
{
  timestamp = format(Sys.time(),"%Y%m%d_%H%M%S")
  results_folder = get_package_folder(paste("results", timestamp, sep="_"))
  dir.create(results_folder)

  results_folder
}
