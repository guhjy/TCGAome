###
## Dataset descriptive analysis
###
descriptive.analysis <- function(X, Y, Z){

  flog.info("Running descriptive analysis...")

  output_dir = paste(RESULTS_FOLDER, "descriptive_analysis", sep="/")
  dir.create(output_dir)

  plot_samples_barplot(Z, output_dir=output_dir, file_name="samples.barplot.png")

  flog.info("Descriptive analysis finished.")
}


###
## Principal Component Analysis
###
pca.analysis <- function(X, Y, Z){

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
hclust.analysis <- function(X, Y, Z){

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
rcaa.analysis <- function(X, Y, Z)
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
rgcca.analysis <- function (X, Y, Z)
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
spls.analysis <- function(X, Y, Z, topN=10)
{

  flog.info("Running sparse Partial Least Squares (sPLS) analysis...")

  output_dir = paste(RESULTS_FOLDER, "spls", sep="/")
  dir.create(output_dir)

  # adds a suffix to gene names to differentiate the dataset of origin, we'll have to remove for posterior analysis
  colnames(X) = paste(colnames(X), ".X", sep="")
  colnames(Y) = paste(colnames(Y), ".Y", sep="")

  # Runs sPLS
  spls_result <- spls(as.matrix(X), as.matrix(Y), ncomp = 3, mode = 'regression')

  # plots the individuals
  mixomics_plot_individuals(data=spls_result, names=Z$Tumor_type, output_dir=output_dir, file_name="samples.png")
  # plots the variables
  mixomics_plot_variables(spls_result, output_dir=output_dir, file_name="variables.png")
  # plots variables in a network
  mixomics_plot_network(spls_result, output_dir=output_dir, file_name="network.png")
  # plots heatmap
  mixomics_plot_heatmap(spls_result, output_dir=output_dir, file_name="heatmap.png")

  # Gets the coordinates of all variables in the correlation plot and calculates the distance to the origin
  X_coords = as.data.frame(cor(spls_result$X,
                    spls_result$variates$X + spls_result$variates$Y,
                    use = "pairwise"))
  names(X_coords) = c("coord1", "coord2", "coord3")
  X_coords$variable = gsub("\\.X", "", row.names(X_coords))
  X_coords$distance2origin = ((X_coords$coord1 ^ 2) + (X_coords$coord2 ^ 2) +(X_coords$coord3 ^ 2))**(1/2)
  Y_coords = as.data.frame(cor(spls_result$Y,
                    spls_result$variates$X + spls_result$variates$Y,
                    use = "pairwise"))
  names(Y_coords) = c("coord1", "coord2", "coord3")
  Y_coords$variable = gsub("\\.Y", "", colnames(spls_result$Y))
  Y_coords$distance2origin = ((Y_coords$coord1 ^ 2) + (Y_coords$coord2 ^ 2) +(Y_coords$coord3 ^ 2))**(1/2)

  # Stores the data
  write.table(X_coords, file = paste(output_dir, "correlation_plot_X_variables_coords.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  write.table(Y_coords, file = paste(output_dir, "correlation_plot_Y_variables_coords.txt", sep="/"), quote=F, row.names=F, sep = "\t")

  # Selects top 20 variables with highest correlations (distance to the origin in the correlaion plot)
  X_selected_variables = X_coords[order(X_coords$distance2origin, decreasing = T)[1:20], ]$variable
  Y_selected_variables = Y_coords[order(Y_coords$distance2origin, decreasing = T)[1:20], ]$variable
  selected_variables = union(X_selected, Y_selected)  # avoids genes identified in both datasets

  # Writes selected variables
  write.table(X_selected_variables, file = paste(output_dir, "RNAseq_selected_variables.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  write.table(Y_selected_variables, file = paste(output_dir, "RPPA_selected_variables.txt", sep="/"), quote=F, row.names=F, sep = "\t")
  write.table(selected_variables, file = paste(output_dir, "selected_variables.txt", sep="/"), quote=F, row.names=F, sep = "\t")

  flog.info("sPLS finished.")

  # return
  list(results=spls_result, selected_variables=selected_variables, selected_genes=X_selected_variables, selected_proteins=Y_selected_variables)
}

###
### Co-Inertia Analysis
###
mcia.analysis <- function(X, Y, Z, topN=5, cia.nf=5)
{

  flog.info("Running Multiple Co-Inertia Analysis (MCIA) ...")

  output_dir = paste(RESULTS_FOLDER, "mcia", sep="/")
  dir.create(output_dir)

  #col.histological.type = as.numeric(Z$histological_type)
  #col.histological.type[col.histological.type == 1] <- red500
  #col.histological.type[col.histological.type == 2] <- green500
  #col.histological.type[col.histological.type == 3] <- indigo500
  #col.histological.type[col.histological.type == 4] <- amber500
  #col.histological.type[col.histological.type == 5] <- yellow500
  #col.histological.type[col.histological.type == 6] <- purple500

  #col.disease.stage = as.numeric(Z$neoplasm_diseasestage)
  #col.disease.stage[col.disease.stage == 1] <- red500
  #col.disease.stage[col.disease.stage == 2] <- green500
  #col.disease.stage[col.disease.stage == 3] <- indigo500
  #col.disease.stage[col.disease.stage == 4] <- amber500
  #col.disease.stage[col.disease.stage == 5] <- yellow500
  #col.disease.stage[col.disease.stage == 6] <- purple500
  #col.disease.stage[col.disease.stage == 7] <- brown500
  #col.disease.stage[col.disease.stage == 8] <- cyan500
  #col.disease.stage[col.disease.stage == 9] <- lime500
  #col.disease.stage[col.disease.stage == 10] <- black
  #col.disease.stage[is.na(col.disease.stage)] <- grey500

  matrices_list = list(as.matrix(t(X)), as.matrix(t(Y)))
  sapply(matrices_list, dim)
  all(apply((x <- sapply(matrices_list, colnames)), 2, function(y)identical(y, x[,1])))

  # Multiple Co-Inertia Analysis
  mcia_result <- mcia(matrices_list, nsc=F, cia.nf = cia.nf)

  # PLot all results
  mcia_plot(mcia_result, phenotype=Z$Tumor_type, output_dir=output_dir, file_name="visualizations.png")

  # Selects top 5 positive and negative associations and plots them, that is 20 variables
  mcia_selected_variables = topVar(mcia_result, topN = topN)
  mcia_selected_genes = gsub("\\.df1", "", as.character(unlist(mcia_selected_variables[, 1:2])))
  mcia_selected_proteins = gsub("\\.df2", "", as.character(unlist(mcia_selected_variables[, 3:4])))
  write.table(mcia_selected_variables, paste(output_dir, "selected_variables.txt", sep="/"), quote = FALSE, sep = "\t", row.names = FALSE)
  mcia_selected_variables = as.character(unlist(mcia_selected_variables))

  # Plots selected variables
  mcia_plot_variables(mcia_result, mcia_selected_variables, output_dir=output_dir, file_name="topN.variables.png")

  # Removes the suffix that identifies a variable as gene or protein
  mcia_selected_variables = gsub("\\.df.", "", mcia_selected_variables)

  flog.info("MCIA finished.")

  list(results=mcia_result, selected_variables=mcia_selected_variables, selected_genes=mcia_selected_genes, selected_proteins=mcia_selected_proteins)
}
