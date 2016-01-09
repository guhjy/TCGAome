#library(treemap)
#library(gridExtra)
#library(ggplot2)
#library(scales)
#library(grid)
#library(VennDiagram)
#library(ggbiplot)
#library(gplots)
#library(Cairo)

# Plots enriched GO terms in a scatter plot after clustering and filtering
plot_scatter <-function(cluster_representatives, output_dir)
{
  # Plot GO terms
  go_viz_names <- c("term_ID", "name", "plot_X","plot_Y","relative_size","qvalue", "annotated_genes", "found_genes", "expected_genes")
  go_viz_data <- cbind(
    cluster_representatives$GO,
    cluster_representatives$name,
    cluster_representatives$x,
    cluster_representatives$y,
    cluster_representatives$size,
    cluster_representatives$qvalue,
    cluster_representatives$annotated_genes,
    cluster_representatives$found_genes,
    cluster_representatives$expected_genes
  )
  go_viz_data = as.data.frame(go_viz_data)
  names(go_viz_data) = go_viz_names
  go_viz_data$plot_X <- as.numeric( as.character(go_viz_data$plot_X) );
  go_viz_data$plot_Y <- as.numeric( as.character(go_viz_data$plot_Y) );
  go_viz_data$relative_size <- as.numeric( as.character(go_viz_data$relative_size) );
  go_viz_data$qvalue <- as.numeric( as.character(go_viz_data$qvalue) );

  Cairo(file= paste(output_dir, "GO_scatterplot.png", sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    p1 <- ggplot( data = go_viz_data );
    p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = qvalue, size = relative_size), alpha = I(0.6) )
    p1 <- p1 + scale_colour_gradientn( "Significance",  colours = c("red", "green"))
    p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = relative_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) ))
    p1 <- p1 + scale_size("Frequency", range=c(5, 20), breaks = c(0.00, 0.02, 0.04), labels = c("0%", "2%", "4%"))
    p1 <- p1 + geom_text(
      data = go_viz_data [ go_viz_data$qvalue < 0.01, ],
      aes(
        plot_X,
        plot_Y,
        label = paste(paste(term_ID, paste(round(relative_size*100, 2), "%"), sep=", "), paste(expected_genes, found_genes, sep=" / "), sep=", ")),
      colour = I(alpha("black", 0.85)),
      size = 3,
      nudge_y = -0.04, check_overlap = F);
    p1 <- p1 + labs (y = "Y", x = "X");
    p1 <- p1 + theme(legend.key = element_blank()) ;
    one.x_range = max(go_viz_data$plot_X) - min(go_viz_data$plot_X);
    one.y_range = max(go_viz_data$plot_Y) - min(go_viz_data$plot_Y);
    p1 <- p1 + xlim(min(go_viz_data$plot_X)-one.x_range/10,max(go_viz_data$plot_X)+one.x_range/10);
    p1 <- p1 + ylim(min(go_viz_data$plot_Y)-one.y_range/10,max(go_viz_data$plot_Y)+one.y_range/10);

  })
  dev.off()

}

# Plots elements in a table
plot_table <- function(cluster_representatives, output_dir){

  cluster_representatives_aux = cluster_representatives[ , c("GO", "name", "qvalue", "size", "expected_genes", "found_genes")]
  cluster_representatives_aux$size = paste(round(cluster_representatives_aux$size*100, 2), "%")
  cluster_representatives_aux$qvalue = round(cluster_representatives_aux$qvalue, 4)

  base_size = 9
  ttheme = ttheme_default(
    core = list(fg_params = list(parse = FALSE,
                                 col = "black",
                                 fontsize=base_size)),
    colhead = list(fg_params = list(parse = FALSE,
                                    fontface = 2L,
                                    fontsize=base_size)),
    rowhead = list(fg_params = list(parse = FALSE,
                                    fontface = 3L,
                                    fontsize=base_size)))

  Cairo(file= paste(output_dir, "GO_terms.png", sep = "/"), type="png", units="in", width=10, height=10, pointsize=12, dpi=600)
  print ({
    tbl <- tableGrob(cluster_representatives_aux, rows=NULL, theme = ttheme)
    grid.draw(tbl)
  })
  dev.off()
}

# Plots GO term clusters in a treemap
plot_treemap <-function(enrichment_results, output_dir)
{

  go_viz_names <- c("term_ID","description","freqInDbPercent","representative");
  go_viz_data <- cbind(
    enrichment_results$GO,
    enrichment_results$name,
    enrichment_results$size,
    enrichment_results$cluster_name
  )

  stuff <- data.frame(go_viz_data);
  names(stuff) <- go_viz_names;

  stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );

  # check the tmPlot command documentation for all possible parameters - there are a lot more
  Cairo(file= paste(output_dir, "GO_treemap.png", sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  treemap(
    stuff,
    index = c("representative","description"),
    vSize = "freqInDbPercent",
    type = "categorical",
    vColor = "representative",
    title = "Gene Ontology treemap",
    inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
    lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
    bg.labels = "#CCCCCCAA",     # define background color of group labels
    # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
    position.legend = "none",
    align.labels=list(c("left", "top"), c("center", "center"))
  )
  dev.off()
}

# Plots GO graph
plot_graph <- function (cluster_representatives, GOdata, output_dir)
{
  scores = as.numeric(cluster_representatives$pvalue)
  names(scores) = cluster_representatives$GO
  Cairo(file=paste(output_dir, "GO_graph.png", sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  showSigOfNodes(GOdata,
                 scores,
                 firstSigNodes = min(10, length(cluster_representatives$GO)),
                 wantedNodes = cluster_representatives$GO,
                 useInfo = 'all')
  dev.off()
}


# Plots Venn diagram
plot_triple_venn <- function(set1, set2, set3, names, output_dir, file_name){
  source("go_plots.r")
  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    draw.triple.venn(area1 = length(set1),
                     area2 = length(set2),
                     area3 = length(set3),
                     n12 = length(intersect(set1, set2)),
                     n23 = length(intersect(set2, set3)),
                     n13 = length(intersect(set1, set3)),
                     n123 = length(intersect(set1, intersect(set2, set3))),
                     category = names, lty = "blank",
                     fill = c(red500, green500, indigo500))
  })
  dev.off()
}

# Plot simple scatter
plot_simple_scatter <- function(matrix, groups, output_dir, file_name){
  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    ggbiplot(matrix, obs.scale = 1, var.scale = 1,
             groups = groups, ellipse = TRUE,
             circle = TRUE, var.axes = FALSE, color = "blue") +
      scale_color_discrete(name = '') +
      theme(legend.direction = 'horizontal', legend.position = 'top')
  })
  dev.off()
}

# Plots simple heatmap
plot_heatmap <- function(matrix, colors, output_dir, file_name){

  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    heatmap.2(matrix, col=redgreen(75), scale="column", key=T, keysize=1.5,
              density.info="none", trace="none",cexCol=0.9, labRow=NA, RowSideColors=colors)
  })
  dev.off()

}

# Plots individuals
mixomics_plot_individuals <- function (data, names, colors, output_dir, file_name){

  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    plotIndiv(data, comp = 1:2, ind.names = names, col = colors)
  })
  dev.off()
}

#Plots variables
mixomics_plot_variables <- function(data, output_dir, file_name){

  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    mixOmics::plotVar(data, comp = 1:2, cutoff = 0.5, Y.label = "Comp 2", X.label = "Comp 1",
                      cex = c(0.8, 0.8))
  })
  dev.off()
}

# Plots a network of features
mixomics_plot_network <- function(spls_result, output_dir, file_name){

  ## By setting keep.var = TRUE, we only display the variables selected by sPLS
  ## on dimensions 1 and 2
  color.edge <- colorRampPalette(c("darkgreen", "green", "yellow", "red", "darkred"))
  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    network(spls_result, comp = 1:2, keep.var = FALSE, shape.node = c("rectangle", "rectangle"),
            color.node = c("white", "pink"), color.edge = color.edge(10), alpha = 3, threshold=0.6)
  })
  dev.off()
}

# Plots a heatmap
mixomics_plot_heatmap <- function (spls_result, output_dir, file_name){

  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    cim(spls_result, comp = 1:3, xlab = "proteins", ylab = "genes", margins = c(5, 6))
  })
  dev.off()
}

# Main MCIA composite plot
mcia_plot <- function(mcia_result, phenotype, output_dir, file_name){

  col.variables = c(indigo500, amber500)
  # TODO: create a color palette with as many colors as phenotypes received
  col.tumor.type = phenotype
  col.tumor.type[col.tumor.type == "BRCA"] <- red500
  col.tumor.type[col.tumor.type == "OV"] <- green500

  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    plot(mcia_result, axes=1:2, phenovec=phenotype, sample.lab=FALSE, df.color=col.variables, sample.color=col.tumor.type,
         sample.legend=FALSE, gene.nlab = 5, df.pch=c(21,22))
  })
  dev.off()
}

mcia_plot_variables <- function(mcia_result, mcia_selected_variables, output_dir, file_name){

  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    omicade4::plotVar(mcia_result, mcia_selected_variables, var.col=red500, var.lab=TRUE, bg.var.col="grey")
  })
  dev.off()
}

plot_samples_barplot <- function(Z, output_dir, file_name){

  Cairo(file=paste(output_dir, file_name, sep = "/"), type="png", units="in", width=10, height=7, pointsize=12, dpi=600)
  print({
    q <- qplot(x=Z$Tumor_type, data=Z, geom="bar", fill = Z$vital_status, ylab = "Count", xlab="Tumor type")
    q  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text( stat='count', aes(label=..count..), vjust=-1, size = 3, hjust = 0.5, position="stack")
  })
  dev.off()
}
