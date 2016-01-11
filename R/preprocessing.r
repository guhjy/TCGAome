
###
## Data preprocessing
###
preprocess_data <- function(X, Y, correlation.thr = 0.7)
{
  flog.info("Starts preprocessing...")

  # Prepares matrices for mixOmics
  X = as.data.frame(X)
  Y = as.data.frame(Y)

  flog.info("No. of RNAseq variables before preprocessing: %s", dim(X)[2])
  flog.info("No. of RPPA variables before preprocessing: %s", dim(Y)[2])

  # Scales matrices
  X = scale(X)
  Y = scale(Y)

  # Remove variables with zero-variance (20169 genes and 142 proteins)
  #variance.thr = 0.01
  X = X[, apply(X, 2, function(x) 0 != var(if (is.factor(x)) as.integer(x) else x)), drop=F]
  Y = Y[, apply(Y, 2, function(x) 0 != var(if (is.factor(x)) as.integer(x) else x)), drop=F]
  X = X[, -nearZeroVar(X)$Position]
  #Y = Y[, -nearZeroVar(Y)$Position]  # there are no nearzero variables in Y

  flog.info("No. of RNAseq variables after removing zero and near zero variance variables: %s", dim(X)[2])
  flog.info("No. of RPPA variables after removing zero and near zero variance variables: %s", dim(Y)[2])

  # Remove variables having any missing value
  # TODO: evaluate other approaches to deal with missing values
  X = X[, apply(X, 2, function(x)!any(is.na(x))), drop=F]
  Y = Y[, apply(Y, 2, function(x)!any(is.na(x))), drop=F]

  flog.info("No. of RNAseq variables after removing variables with missing values: %s", dim(X)[2])
  flog.info("No. of RPPA variables after removing variables with missing values: %s", dim(Y)[2])

  # Removes correlated variables using a threshold of 0.7 (14179 genes and 121 proteins)
  #correlation.thr = 0.7  ## This value is in the literature from Le Cao
  X.cor = cor(X)
  X.cor[!lower.tri(X.cor)] <- 0 # sets upper triangle and diagonal to 0
  X <- X[,!apply(X.cor,2,function(x) any(x > correlation.thr))]
  Y.cor = cor(Y)
  Y.cor[!lower.tri(Y.cor)] <- 0 # sets upper triangle and diagonal to 0
  Y <- Y[,!apply(Y.cor,2,function(x) any(x > correlation.thr))]

  flog.info("Correlation threshold: %s", correlation.thr)
  flog.info("No. of RNAseq variables after removing duplicate variables with high correlation: %s", dim(X)[2])
  flog.info("No. of RPPA variables after removing duplicate variables with high correlation: %s", dim(Y)[2])

  flog.info("Preprocessing finished.")

  # returns preprocessed matrices X, Y
  list(X=X, Y=Y)
}
