###
## Dataset download and preparation
###
get_tcga_data <- function (tumor_types)
{
  flog.info("Starting download for tumor types: %s", paste(tumor_types, collapse=", "))
  firehose_datasets = getFirehoseDatasets()
  #firehose_dates= getFirehoseRunningDates()
  FIREHOSE_DATE = "20150402"  # we fix it so data is consistent during development
  FIREHOSE_DIR <<- get_package_folder("inst/firehose")

  datasets = list()     # whole datasets
  matrices = list()
  X = NULL
  Y = NULL
  Z = NULL

  #tumor_types = c("BRCA", "OV")

  for (tumor_type in tumor_types)
  {
    flog.info("Download data for %s", tumor_type)
    # downloads the data
    datasets[[tumor_type]] = getFirehoseData (dataset=tumor_type, runDate=FIREHOSE_DATE, RNAseq2_Gene_Norm = TRUE, RPPA=TRUE, destdir=FIREHOSE_DIR)

    flog.info("Loading data for %s", tumor_type)
    # gets data matrices
    matrices[[tumor_type]]$X = datasets[[tumor_type]]@RNASeq2GeneNorm
    matrices[[tumor_type]]$Y = datasets[[tumor_type]]@RPPAArray[[1]]@DataMatrix
    matrices[[tumor_type]]$Z = as.data.frame(datasets[[tumor_type]]@Clinical)
    if (grepl("-.-(V|C|NA)$",rownames(matrices[[tumor_type]]$Y))){
      rownames(matrices[[tumor_type]]$Y) = gsub("-.-(V|C|NA)$", "", rownames(matrices[[tumor_type]]$Y))  # antobodies names include a suffix for OV dataset...
    }

    # Normalize sample ids
    matrices[[tumor_type]] = normalize_sample_identifiers(matrices[[tumor_type]]$X,
                                            matrices[[tumor_type]]$Y,
                                            matrices[[tumor_type]]$Z)

    flog.info("No. of %s cancer samples: %s", tumor_type, dim(matrices[[tumor_type]]$X)[2])
    flog.info("No. of %s cancer features for RNAseq: %s", tumor_type, dim(matrices[[tumor_type]]$X)[1])
    flog.info("No. of %s cancer features for RPPA: %s", tumor_type, dim(matrices[[tumor_type]]$Y)[1])
    flog.info("No. of %s cancer features for clinical data: %s", tumor_type, dim(matrices[[tumor_type]]$Z)[2])
    #matrices[[tumor_type]]$Z = matrices[[tumor_type]]$Z

    matrices[[tumor_type]] = get_RPPA_annotations(tumor_type, matrices[[tumor_type]])

    # adds the tumor type in the clinical data matrix
    matrices[[tumor_type]]$Z = cbind(matrices[[tumor_type]]$Z, Tumor_type=c(tumor_type))

    # unlists matrices of each data type for later processing
    if (is.null(X)){
      X = matrices[[tumor_type]]$X
    } else {
      shared_features = intersect(rownames(X), rownames(matrices[[tumor_type]]$X))
      X = cbind2(X[rownames(X) %in% shared_features, ],
                 matrices[[tumor_type]]$X[rownames(matrices[[tumor_type]]$X) %in% shared_features, ])
    }
    if (is.null(Y)){
      Y = matrices[[tumor_type]]$Y
    } else {
      shared_features = intersect(rownames(Y), rownames(matrices[[tumor_type]]$Y))
      Y = cbind2(Y[rownames(Y) %in% shared_features, ],
                 matrices[[tumor_type]]$Y[rownames(matrices[[tumor_type]]$Y) %in% shared_features, ])
    }
    if (is.null(Z)){
      Z = matrices[[tumor_type]]$Z
    } else {
      Z = subset(Z, select=intersect(colnames(Z), colnames(matrices[[tumor_type]]$Z)))
      matrices[[tumor_type]]$Z = subset(matrices[[tumor_type]]$Z, select=intersect(colnames(Z), colnames(matrices[[tumor_type]]$Z)))
      Z = rbind2(Z, matrices[[tumor_type]]$Z)
    }
  }

  # 18147 gene expression variables + 133 protein expression variables for 407 samples
  flog.info("Total number of features for RNAseq: %s", dim(X)[1])
  flog.info("Total number of features for RPPA: %s", dim(Y)[1])
  flog.info("Total number of features for clinical data: %s", dim(Z)[2])
  flog.info("Total number of samples: %s", dim(X)[2])

  # Transposes
  X = t(X)
  Y = t(Y)

  # Sorts matrices by sample
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  Z <- Z[order(rownames(Z)),]

  # Writes input matrices
  write.table(X, file=paste(RESULTS_FOLDER, "X_rnaseqnorm.txt", sep="/"), quote = F, sep = "\t", na = "")
  write.table(Y, file=paste(RESULTS_FOLDER, "Y_rppa.txt", sep="/"), quote = F, sep = "\t", na = "")
  write.table(Z, file=paste(RESULTS_FOLDER, "Z_clinical.txt", sep="/"), quote = F, sep = "\t", na = "")

  flog.info("Finished downloading TCGA data.")

  # returns matrices
  list(X=X, Y=Y, Z=Z)
}

RPPA_ANNOTATIONS = list(
  "ACC" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/ACC/20150821/gdac.broadinstitute.org_ACC.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "BLCA" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/BLCA/20150821/gdac.broadinstitute.org_BLCA.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "BRCA" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/BRCA/20150821/gdac.broadinstitute.org_BRCA.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "CESC" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/CESC/20150821/gdac.broadinstitute.org_CESC.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "CHOL" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/CHOL/20150821/gdac.broadinstitute.org_CHOL.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "COAD" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/COAD/20150821/gdac.broadinstitute.org_COAD.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "COADREAD" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/COADREAD/20150821/gdac.broadinstitute.org_COADREAD.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "DLBC" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/DLBC/20150821/gdac.broadinstitute.org_DLBC.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "ESCA" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/ESCA/20150821/gdac.broadinstitute.org_ESCA.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "GBM" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/GBM/20150821/gdac.broadinstitute.org_GBM.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "GBMLGG" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/GBMLGG/20150821/gdac.broadinstitute.org_GBMLGG.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "HNSC" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/HNSC/20150821/gdac.broadinstitute.org_HNSC.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "KICH" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/KICH/20150821/gdac.broadinstitute.org_KICH.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "KIPAN" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/KIPAN/20150821/gdac.broadinstitute.org_KIPAN.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "KIRC" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/KIRC/20150821/gdac.broadinstitute.org_KIRC.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "KIRP" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/KIRP/20150821/gdac.broadinstitute.org_KIRP.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "LGG" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/LGG/20150821/gdac.broadinstitute.org_LGG.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "LIHC" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/LIHC/20150821/gdac.broadinstitute.org_LIHC.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "LUAD" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/LUAD/20150821/gdac.broadinstitute.org_LUAD.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "LUSC" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/LUSC/20150821/gdac.broadinstitute.org_LUSC.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "MESO" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/MESO/20150821/gdac.broadinstitute.org_MESO.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "OV"="http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/OV/20150821/gdac.broadinstitute.org_OV.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "PAAD" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/PAAD/20150821/gdac.broadinstitute.org_PAAD.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "PCPG" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/PCPG/20150821/gdac.broadinstitute.org_PCPG.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "PRAD" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/PRAD/20150821/gdac.broadinstitute.org_PRAD.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "READ" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/READ/20150821/gdac.broadinstitute.org_READ.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "SARC" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/SARC/20150821/gdac.broadinstitute.org_SARC.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "SKCM" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/SKCM/20150821/gdac.broadinstitute.org_SKCM.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "STAD" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/STAD/20150821/gdac.broadinstitute.org_STAD.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "STES" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/STES/20150821/gdac.broadinstitute.org_STES.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "TGCT" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/TGCT/20150821/gdac.broadinstitute.org_TGCT.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "THCA" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/THCA/20150821/gdac.broadinstitute.org_THCA.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "THYM" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/THYM/20150821/gdac.broadinstitute.org_THYM.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "UCEC" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/UCEC/20150821/gdac.broadinstitute.org_UCEC.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz",
  "UCS" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/UCS/20150821/gdac.broadinstitute.org_UCS.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz",
  "UVM" = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/UVM/20150821/gdac.broadinstitute.org_UVM.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz"
  )

# Downloads annotations of RPPA antibodies for a given tumor type
get_RPPA_annotations <- function(tumor_type, matrices){

  flog.info("Loading RPPA annotations for %s", tumor_type)

  #RPPA.annotations.URL.template = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/%TUMOR_TYPE%/20150821/gdac.broadinstitute.org_%TUMOR_TYPE%.RPPA_AnnotateWithGene.Level_3.2015082100.1.0.tar.gz"
  #RPPA.annotations.URL.template_2 = "http://gdac.broadinstitute.org/runs/stddata__2015_08_21/data/%TUMOR_TYPE%/20150821/gdac.broadinstitute.org_%TUMOR_TYPE%.RPPA_AnnotateWithGene.Level_3.2015082100.0.0.tar.gz"
  RPPA.annotations.file.template = paste(FIREHOSE_DIR, "%TUMOR_FOLDER%/%TUMOR_TYPE%.antibody_annotation.txt", sep="/")

  # Downloads the annotations if necessary
  RPPA.annotations.URL = RPPA_ANNOTATIONS[[tumor_type]]
  RPPA.annotations.file = paste(FIREHOSE_DIR, basename(RPPA.annotations.URL), sep="/")
  if (!file.exists(RPPA.annotations.file)){
    download.file(RPPA.annotations.URL, RPPA.annotations.file)
  }

  # Uncompressing annotations
  if (file.exists(RPPA.annotations.file)){
    untar(RPPA.annotations.file, exdir = FIREHOSE_DIR)
    output_matrices = normalize_variable_names(matrices$X, matrices$Y, matrices$Z,
                                               gsub("%TUMOR_FOLDER%",
                                                    gsub(".tar.gz", "",
                                                         basename(RPPA.annotations.file)),
                                                    gsub("%TUMOR_TYPE%",
                                                         tumor_type, RPPA.annotations.file.template)))
  } else {
    # Annotations not found where expected
    flog.error("Cannot find RPPA annotations for tumor type %s.", tumor_type)
  }

  flog.info("Finished loading RPPA annotations.")

  output_matrices
}

# Normalize sample identifiers
normalize_sample_identifiers <- function (X, Y, Z)
{
  flog.info("Normalizing sample identifiers...")

  # gets samples
  gene_expression_samples = colnames(X)
  protein_expression_samples = colnames(Y)

  # Normalizes sample identifiers for X and Y
  # See https://wiki.nci.nih.gov/display/TCGA/TCGA+Barcode for details on sample identifiers
  # Check https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=Sample%20type to identify tumor and normal samples
  # removes the part in the sample id related with the data type
  protein_expression_samples_processed = unlist(lapply(protein_expression_samples, substr, start=1, stop=16))
  colnames(Y) = protein_expression_samples_processed
  gene_expression_samples_processed = unlist(lapply(gene_expression_samples, substr, start=1, stop=16))
  colnames(X) = gene_expression_samples_processed

  # Normalizes sample identifiers for clinical data (identifiers are of the type "TCGA-MS-A51U", they don't contain the part related with the tissue)
  row.names(Z) = gsub("\\.", "-", toupper(row.names(Z)))

  # Gets common samples (407 samples having gene expression and protein expression)
  common_samples = intersect(protein_expression_samples_processed, gene_expression_samples_processed)

  # Get samples subset based on histology type


  # Removes samples not being common from both data types
  X = subset(X, select= common_samples)
  Y = subset(Y, select= common_samples)

  # TODO: Remove repetitions for the same individual (there are none in this dataset as we have no controls for protein expression)
  # common_samples = uniq(substr(common_samples, 1, 12))

  # Normalizes sample identifiers with clinical data
  # Gets only the part of the sample identifier corresponding to individual so we can match to clinical data ("TCGA-MS-A51U")
  colnames(X)  = substr(colnames(X), 1, 12)
  colnames(Y)  = substr(colnames(Y), 1, 12)
  Z_t = subset(t(Z), select= intersect(colnames(Y), rownames(Z)))
  Z = as.data.frame(t(Z_t))
  rownames(Z)  = substr(rownames(Z), 1, 12)

  #return
  list(X=X, Y=Y, Z=Z)
}

# Normalizes variables into HGNC symbols
normalize_variable_names <- function(X, Y, Z, antibody.annotation)
{

  flog.info("Normalizing feature names...")

  # Normalizes gene identifiers using Ensembl BioMart for Homo sapiens
  # Removes genes from X not having a correct HUGO identifier
  hgnc_symbols = query_biomart(attributes=c("hgnc_symbol"), filters=c("hgnc_symbol"), values=rownames(X))[, 1]
  X = X[ rownames(X) %in% hgnc_symbols, ]

  # Normalizes RPPA antibody identifiers with Firehose annotation to gene (we use BRCA annotations as they are subset of OV)
  # Some antibodies map to several genes -> we will keep only the first in the list
  # Some genes are mapped by several antibodies -> we will keep only the one with the highest average expression
  antibody.annotation = read.csv(antibody.annotation, sep="\t")
  # remove more than one gene to the same antibody
  antibody.annotation$Gene.Name.processed =  gsub(" .*", "", antibody.annotation$Gene.Name)
  antibodies = intersect(rownames(Y), antibody.annotation$Composite.Element.REF)
  Y = Y[ rownames(Y) %in% antibodies, ]

  # remove more than one antibody to the same gene, we keep that one with the greatest average expression
  antibody.annotation = antibody.annotation[antibody.annotation$Composite.Element.REF %in% antibodies, ]
  antibody.annotation$max.avg.exp = apply(Y[rownames(Y) %in% antibodies, ], MARGIN = 1, FUN = function(x){abs(mean(x))})
  antibody.annotation <- antibody.annotation[order(antibody.annotation$Gene.Name.processed, antibody.annotation$max.avg.exp, decreasing=TRUE),]
  antibody.annotation <- antibody.annotation[!duplicated(antibody.annotation$Gene.Name.processed),]
  antibodies = intersect(rownames(Y), antibody.annotation$Composite.Element.REF)
  Y = Y[ rownames(Y) %in% antibodies, ]

  # Replace antibody identifiers with gene symbols
  antibody.annotation <- antibody.annotation[order(antibody.annotation$Composite.Element.REF),]
  Y <- Y[order(rownames(Y)),]
  rownames(Y) = antibody.annotation$Gene.Name.processed

  # return
  list(X=X, Y=Y, Z=Z)
}
