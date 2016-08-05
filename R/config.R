

# Configures loggers to be used
configure_logging <- function(folder) {
    futile.logger::flog.threshold(futile.logger::DEBUG)
    futile.logger::flog.appender(futile.logger::appender.tee(paste(folder, "TCGAome.log", sep = "/")))
}

# Gets the absolute path for the package folder
get_package_folder <- function(relative_path) {
    folder <- file.path(gsub("\\\\", "/", base::system.file(package = "TCGAome")), relative_path)
    dir.create(folder, showWarnings = FALSE, recursive = T)
    futile.logger::flog.debug("Package folder: %s", folder)
    return(folder)
}

# Gets the results folder with a unique timestamp
get_results_folder <- function() {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    results_folder <- get_package_folder(paste("results", timestamp, sep = "_"))
    dir.create(results_folder, showWarnings = FALSE, recursive = T)
    futile.logger::flog.debug("Results folder: %s", results_folder)
    results_folder
}
