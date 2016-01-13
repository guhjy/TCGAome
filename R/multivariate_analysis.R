###
## Dataset descriptive analysis
###
descriptive_analysis <- function(X, Y, Z){

  flog.info("Running descriptive analysis...")

  output_dir = paste(RESULTS_FOLDER, "descriptive_analysis", sep="/")
  dir.create(output_dir)

  plot_samples_barplot(Z, output_dir=output_dir, file_name="samples.barplot.png")

  flog.info("Descriptive analysis finished.")
}


###
## Principal Component Analysis
###
pca_analysis <- function(X, Y, Z){

  flog.info("Running PCA analysis...")

  output_dir = paste(RESULTS_FOLDER, "pca", sep="/")
  dir.create(output_dir)

  # PCA on the protein expression
  Y.pca <- prcomp(Y, center = TRUE, scale. = TRUE)
  plot_simple_scatter(matrix=Y.pca, groups=Z$Tumor_type, output_dir=output_dir, file_name="PCA.Y.png")

  # PCA on the gene expression
  X.pca <- prcomp(X,
                  center = TRUE,
                  scale. = TRUE)
  plot_simple_scatter(matrix=X.pca, groups=Z$Tumor_type, output_dir=output_dir, file_name="PCA.X.png")

  # Joint PCA
  colnames(X) = paste(colnames(X), ".X")
  colnames(Y) = paste(colnames(Y), ".Y")
  XY.pca <- prcomp(cbind2(X,Y),
                   center = TRUE,
                   scale. = TRUE)
  plot_simple_scatter(matrix=XY.pca, groups=Z$Tumor_type, output_dir=output_dir, file_name="PCA.XY.png")

  flog.info("PCA analysis finished.")
}

###
## Hierarchical clustering
###
hclust_analysis <- function(X, Y, Z){

  flog.info("Running hierarchical clustering analysis...")

  output_dir = paste(RESULTS_FOLDER, "hclust", sep="/")
  dir.create(output_dir)

  colnames(X) = paste(colnames(X), ".X")
  colnames(Y) = paste(colnames(Y), ".Y")

  # Plots heatmaps
  plot_heatmap(matrix=as.matrix(X), tumor_types=Z$Tumor_type, output_dir=output_dir, file_name="hclust.X.png")
  plot_heatmap(matrix=as.matrix(Y), tumor_types=Z$Tumor_type, output_dir=output_dir, file_name="hclust.Y.png")
  plot_heatmap(matrix=as.matrix(cbind2(X,Y)), tumor_types=Z$Tumor_type, output_dir=output_dir, file_name="hclust.XY.png")

  flog.info("Hierarchical clustering finished.")
}

###
### Regularized Canonical COrrelation Analysis (rCCA)
###
rcaa_analysis <- function(X, Y, Z)
{
  flog.info("Running Regularized Canonical Correlation Analysis (rCCA)...")

  output_dir = paste(RESULTS_FOLDER, "rCCA", sep="/")
  dir.create(output_dir)

  # Calculates regularization factors with CV (M-folds and leave-one-out methods raise an error)
  # Error in chol.default(Cxx) :
  # the leading minor of order 407 is not positive definite
  #grid1 <- seq(0, 0.2, length = 51)
  #grid2 <- seq(0.0001, 0.2, length = 51)
  cv_score <- tune.rcc(X, Y)  #, grid1 = grid1, grid2 = grid2, plt = FALSE) #, validation = 'loo')
  # Computes rCCA
  rcca_results <- rcc(X, Y, lambda1=cv_score$opt.lambda1, lambda2=cv_score$opt.lambda2)

  # plots the individuals
  mixomics_plot_individuals(data=rcca_results, names=Z$Tumor_type, output_dir=output_dir, file_name="samples.png")

  # plots the variables
  mixomics_plot_variables(rcca_results, output_dir=output_dir, file_name="variables.png")

  # network of relations between variables
  #network(rcca_results, comp = 1:2, threshold = 0.2)

  # heatmap
  #cim(rcca_results, comp = 1:2, xlab = "proteins", ylab = "genes", margins = c(5,6))

  flog.info("rCCA finished.")

  # return
  rcca_results
}

###
## Regularized Generalized CCA
###
rgcca_analysis <- function (X, Y, Z)
{
  flog.info("Running Regularized Generalized Canonical Correlation Analysis (RGCCA)...")

  rgcca.results = rgcca(list(X, Y), tau = 'optimal', ncomp = 3)

  flog.info("RGCCA finished.")

  # return
  rgcca.results
}


###
### Sparse Partial Least Squares (sPLS)
###
spls_analysis <- function(X, Y, Z, topN=5, selection_method="loadings")
{

  flog.info("Running sparse Partial Least Squares (sPLS) analysis...")

  if (!selection_method %in% c("correlation", "loadings")){
    flog.error("Selection method '%s' non supported.", selection_method)
    stop("Selection method non supported!")
  }

  output_dir = paste(RESULTS_FOLDER, "spls", sep="/")
  dir.create(output_dir)

  # adds a suffix to gene names to differentiate the dataset of origin, we'll have to remove for posterior analysis
  colnames(X) = paste(colnames(X), ".X", sep="")
  colnames(Y) = paste(colnames(Y), ".Y", sep="")

  # Runs sPLS
  spls_result <- spls(as.matrix(X), as.matrix(Y), ncomp = 3, mode = 'regression')

  # plots the individuals
  try(mixomics_plot_individuals(data=spls_result, names=Z$Tumor_type, output_dir=output_dir, file_name="samples.png"), silent = T)
  # plots the variables
  try(mixomics_plot_variables(spls_result, output_dir=output_dir, file_name="variables.png"), silent = T)
  # plots variables in a network
  try(mixomics_plot_network(spls_result, output_dir=output_dir, file_name="network.png"), silent = T)
  # plots heatmap
  try(mixomics_plot_heatmap(spls_result, output_dir=output_dir, file_name="heatmap.png"), silent = T)

  # Gets the coordinates of all variables in the correlation plot and calculates the distance to the origin
  X_coords = as.data.frame(cor(spls_result$X, spls_result$variates$X + spls_result$variates$Y, use = "pairwise"))
  names(X_coords) = c("coord1", "coord2", "coord3")
  X_coords$variable = gsub("\\.X", "", row.names(X_coords))
  X_coords$distance2origin = ((X_coords$coord1 ^ 2) + (X_coords$coord2 ^ 2) +(X_coords$coord3 ^ 2))**(1/2)
  Y_coords = as.data.frame(cor(spls_result$Y, spls_result$variates$X + spls_result$variates$Y, use = "pairwise"))
  names(Y_coords) = c("coord1", "coord2", "coord3")
  Y_coords$variable = gsub("\\.Y", "", colnames(spls_result$Y))
  Y_coords$distance2origin = ((Y_coords$coord1 ^ 2) + (Y_coords$coord2 ^ 2) +(Y_coords$coord3 ^ 2))**(1/2)

  # Stores the data
  write.table(X_coords, file = paste(output_dir, "correlation_plot_X_variables_coords.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  write.table(Y_coords, file = paste(output_dir, "correlation_plot_Y_variables_coords.txt", sep="/"), quote=F, row.names=F, sep = "\t")

  # Reads variables loadings
  variable_loadings_comp1 = mixOmics::selectVar(spls_result, comp=1)
  variable_loadings_comp2 = mixOmics::selectVar(spls_result, comp=2)
  variable_loadings_comp3 = mixOmics::selectVar(spls_result, comp=3)

  # Stores the loadings data
  X_loadings_comp1 = variable_loadings_comp1$X$value
  names(X_loadings_comp1) = c("comp1")
  X_loadings_comp1$variable = row.names(X_loadings_comp1)
  X_loadings_comp2 = variable_loadings_comp2$X$value
  names(X_loadings_comp2) = c("comp2")
  X_loadings_comp2$variable = row.names(X_loadings_comp2)
  X_loadings_comp3 = variable_loadings_comp3$X$value
  names(X_loadings_comp3) = c("comp3")
  X_loadings_comp3$variable = row.names(X_loadings_comp3)
  X_loadings = merge.data.frame(X_loadings_comp1, X_loadings_comp2, by = "variable")
  X_loadings = merge.data.frame(X_loadings, X_loadings_comp3, by = "variable")
  write.table(X_loadings, file = paste(output_dir, "X_variables_loadings.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  Y_loadings_comp1 = variable_loadings_comp1$Y$value
  names(Y_loadings_comp1) = c("comp1")
  Y_loadings_comp1$variable = row.names(Y_loadings_comp1)
  Y_loadings_comp2 = variable_loadings_comp2$Y$value
  names(Y_loadings_comp2) = c("comp2")
  Y_loadings_comp2$variable = row.names(Y_loadings_comp2)
  Y_loadings_comp3 = variable_loadings_comp3$Y$value
  names(Y_loadings_comp3) = c("comp3")
  Y_loadings_comp3$variable = row.names(Y_loadings_comp3)
  Y_loadings = merge.data.frame(Y_loadings_comp1, Y_loadings_comp2, by = "variable")
  Y_loadings = merge.data.frame(Y_loadings, Y_loadings_comp3, by = "variable")
  write.table(Y_loadings, file = paste(output_dir, "Y_variables_loadings.txt", sep="/"), quote=F, row.names=F, sep = "\t")

  if (selection_method == "correlation"){
    # Selects top N variables with highest correlations (distance to the origin in the correlaion plot)
    X_selected_variables = X_coords[order(X_coords$distance2origin, decreasing = T)[1:topN], ]$variable
    Y_selected_variables = Y_coords[order(Y_coords$distance2origin, decreasing = T)[1:topN], ]$variable
  } else if (selection_method == "loadings"){
    # Selects top N variables with highest loadings, maximize variance in samples
    # select topN t each component for positive vaues and for negative values separately
    X_top_negative_loadings_comp1 = gsub(".X", "", variable_loadings_comp1$X$name[variable_loadings_comp1$X$value$value.var <  0][1:topN])
    X_top_positive_loadings_comp1 = gsub(".X", "", variable_loadings_comp1$X$name[variable_loadings_comp1$X$value$value.var >= 0][1:topN])
    Y_top_negative_loadings_comp1 = gsub(".Y", "", variable_loadings_comp1$Y$name[variable_loadings_comp1$Y$value <  0][1:topN])
    Y_top_positive_loadings_comp1 = gsub(".Y", "", variable_loadings_comp1$Y$name[variable_loadings_comp1$Y$value >= 0][1:topN])
    X_top_negative_loadings_comp2 = gsub(".X", "", variable_loadings_comp2$X$name[variable_loadings_comp2$X$value <  0][1:topN])
    X_top_positive_loadings_comp2 = gsub(".X", "", variable_loadings_comp2$X$name[variable_loadings_comp2$X$value >= 0][1:topN])
    Y_top_negative_loadings_comp2 = gsub(".Y", "", variable_loadings_comp2$Y$name[variable_loadings_comp2$Y$value <  0][1:topN])
    Y_top_positive_loadings_comp2 = gsub(".Y", "", variable_loadings_comp2$Y$name[variable_loadings_comp2$Y$value >= 0][1:topN])
    X_top_negative_loadings_comp3 = gsub(".X", "", variable_loadings_comp3$X$name[variable_loadings_comp3$X$value <  0][1:topN])
    X_top_positive_loadings_comp3 = gsub(".X", "", variable_loadings_comp3$X$name[variable_loadings_comp3$X$value >= 0][1:topN])
    Y_top_negative_loadings_comp3 = gsub(".Y", "", variable_loadings_comp3$Y$name[variable_loadings_comp3$Y$value <  0][1:topN])
    Y_top_positive_loadings_comp3 = gsub(".Y", "", variable_loadings_comp3$Y$name[variable_loadings_comp3$Y$value >= 0][1:topN])

    X_top_loadings_comp1 = union(X_top_negative_loadings_comp1, X_top_positive_loadings_comp1)
    Y_top_loadings_comp1 = union(Y_top_negative_loadings_comp1, Y_top_positive_loadings_comp1)
    X_top_loadings_comp2 = union(X_top_negative_loadings_comp2, X_top_positive_loadings_comp2)
    Y_top_loadings_comp2 = union(Y_top_negative_loadings_comp2, Y_top_positive_loadings_comp2)
    X_top_loadings_comp3 = union(X_top_negative_loadings_comp3, X_top_positive_loadings_comp3)
    Y_top_loadings_comp3 = union(Y_top_negative_loadings_comp3, Y_top_positive_loadings_comp3)

    X_selected_variables = union(X_top_loadings_comp1, union(X_top_loadings_comp2, X_top_loadings_comp3))
    Y_selected_variables = union(Y_top_loadings_comp1, union(Y_top_loadings_comp2, Y_top_loadings_comp3))
  }
  selected_variables = union(X_selected_variables, Y_selected_variables)  # avoids genes identified in both datasets

  # Writes selected variables
  write.table(X_selected_variables, file = paste(output_dir, "RNAseq_selected_variables.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  write.table(Y_selected_variables, file = paste(output_dir, "RPPA_selected_variables.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  write.table(selected_variables, file = paste(output_dir, "selected_variables.txt", sep="/"), quote=F, row.names=F, sep = "\t")

  # TODO: plots selected variables

  flog.info("%d variables selected from RNAseq on sPLS results.", length(X_selected_variables))
  flog.info("%d variables selected from RPPA on sPLS results.", length(Y_selected_variables))
  flog.info("%d variables selected on sPLS results.", length(selected_variables))
  flog.info("sPLS finished.")

  # return
  list(results=spls_result, selected_variables=selected_variables, selected_genes=X_selected_variables, selected_proteins=Y_selected_variables)
}

###
### Co-Inertia Analysis
###
mcia_analysis <- function(X, Y, Z, topN=5, cia.nf=5)
{

  flog.info("Running Multiple Co-Inertia Analysis (MCIA) ...")

  output_dir = paste(RESULTS_FOLDER, "mcia", sep="/")
  dir.create(output_dir)

  matrices_list = list(as.matrix(t(X)), as.matrix(t(Y)))
  sapply(matrices_list, dim)
  all(apply((x <- sapply(matrices_list, colnames)), 2, function(y)identical(y, x[,1])))

  # Multiple Co-Inertia Analysis
  mcia_result <- mcia(matrices_list, nsc=F, cia.nf = cia.nf)

  # Extracts Co-Inertia plot coords
  coinertia_variable_coords = mcia_result$mcoa$axis
  coinertia_X_variable_coords = coinertia_variable_coords[grep("\\.df1", row.names(coinertia_variable_coords)), ]
  coinertia_Y_variable_coords = coinertia_variable_coords[grep("\\.df2", row.names(coinertia_variable_coords)), ]
  coinertia_X_variable_coords$variable = gsub("\\.df.", "", row.names(coinertia_X_variable_coords))
  coinertia_Y_variable_coords$variable = gsub("\\.df.", "", row.names(coinertia_Y_variable_coords))
  write.table(coinertia_X_variable_coords, file = paste(output_dir, "coinertia_plot_X_variables_coords.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  write.table(coinertia_Y_variable_coords, file = paste(output_dir, "coinertia_plot_Y_variables_coords.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  coinertia_sample_coords = mcia_result$mcoa$SynVar
  coinertia_sample_coords$sample = row.names(coinertia_sample_coords)
  write.table(coinertia_sample_coords, file = paste(output_dir, "coinertia_sample_coords.txt", sep="/"), quote=F, row.names=F, sep = "\t")

  # PLot all results
  try(mcia_plot(mcia_result, phenotype=Z$Tumor_type, output_dir=output_dir, file_name="visualizations.png"), silent = T)

  # Selects top N positive and negative associations at X and Y and plots them, that is a maximum of 3*2*2*N variables
  selected_variables_1 = topVar(mcia_result, axis=1, topN = topN)
  X_selected_variables_1 = gsub("\\.df1", "", as.character(unlist(selected_variables_1[, 1:2])))
  Y_selected_variables_1 = gsub("\\.df2", "", as.character(unlist(selected_variables_1[, 3:4])))
  selected_variables_1 = as.character(unlist(selected_variables_1[, 1:4]))
  selected_variables_2 = topVar(mcia_result, axis=2, topN = topN)
  X_selected_variables_2 = gsub("\\.df1", "", as.character(unlist(selected_variables_2[, 1:2])))
  Y_selected_variables_2 = gsub("\\.df2", "", as.character(unlist(selected_variables_2[, 3:4])))
  selected_variables_2 = as.character(unlist(selected_variables_2[, 1:4]))
  selected_variables_3 = topVar(mcia_result, axis=3, topN = topN)
  X_selected_variables_3 = gsub("\\.df1", "", as.character(unlist(selected_variables_3[, 1:2])))
  Y_selected_variables_3 = gsub("\\.df2", "", as.character(unlist(selected_variables_3[, 3:4])))
  selected_variables_3 = as.character(unlist(selected_variables_3[, 1:4]))

  selected_variables = union(selected_variables_1, union(selected_variables_2, selected_variables_3))
  X_selected_variables = union(X_selected_variables_1, union(X_selected_variables_2, X_selected_variables_3))
  Y_selected_variables = union(Y_selected_variables_1, union(Y_selected_variables_2, Y_selected_variables_3))

  # Plots selected variables
  try(mcia_plot_variables(mcia_result, selected_variables, output_dir=output_dir, file_name="topN.variables.png"), silent = T)

  # TODO: plot 3rd dimension

  # Removes the suffix that identifies a variable as gene or protein before writing
  selected_variables = gsub("\\.df.", "", selected_variables)

  # Stores selected variables
  write.table(X_selected_variables, file = paste(output_dir, "RNAseq_selected_variables.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  write.table(Y_selected_variables, file = paste(output_dir, "RPPA_selected_variables.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  write.table(selected_variables, file = paste(output_dir, "selected_variables.txt", sep="/"), quote=F, row.names=F, sep = "\t")

  flog.info("%d variables selected from RNAseq on MCIA results.", length(X_selected_variables))
  flog.info("%d variables selected from RPPA on MCIA results.", length(Y_selected_variables))
  flog.info("%d variables selected on MCIA results.", length(selected_variables))
  flog.info("MCIA finished.")

  list(results=mcia_result, selected_variables=selected_variables, selected_genes=X_selected_variables, selected_proteins=Y_selected_variables)
}
