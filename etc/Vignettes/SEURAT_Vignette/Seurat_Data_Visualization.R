#### DATA VISUALIZATION METHODS IN SEURAT ######################################
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
InstallData("pbmc3k")
data("pbmc3k.final")

pbmc3k.final$groups <- sample(c("group1", "group2"),
                              size = ncol(pbmc3k.final),
                              replace = T)
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final
## An object of class Seurat 
## 13714 features across 2638 samples within 1 assay 
## Active assay: RNA (13714 features, 2000 variable features)
##  2 dimensional reductions calculated: pca, umap



#### Visualizations of marker feature expression

## RIDGE PLOTS - from ggridges
    ## Visualize single cell expression distributions in each cluster
    RidgePlot(pbmc3k.final,
              features = features,
              ncol = 2)

## VIOLIN PLOTS
    ## visualize single cell expression distributions in each cluster
    VlnPlot(pbmc3k.final, features = features)
    
    ## Violin Plots can also be split on some variable
    ## Simply add the splitting variable to object metadata
    ## and pass it to the split.by argument
    VlnPlot(pbmc3k.final, 
            features = "percent.mt",
            split.by = "groups")


## FEATURE PLOTS
    ## visualize feature expression in low-dimensional space
    FeaturePlot(pbmc3k.final, features = features)
    
    #### New Additions to FeaturePlot
    ## Plot a legend to map colors to expression levels
    FeaturePlot(pbmc3k.final, 
                features = "MS4A1",
                min.cutoff = 1,
                max.cutoff = 3)


    ## Calculate feature-specific contrast levels based on
    ## quantiles of non-zero expression.
    ## Particularly useful when plotting multiple markers
    FeaturePlot(pbmc3k.final,
                features = c("MS4A1",
                             "PTPRCAP"),
                min.cutoff = "q10",
                max.cutoff = "q90")


    ## Visualize Co-Expression of 2 features simultaneously
    FeaturePlot(pbmc3k.final, features = c("MS4A1",
                                           "CD79A"),
                blend = T)
    
    ## Split visualization to view by groups (replaces FeatureHeatmap)
    FeaturePlot(pbmc3k.final, 
                features = c("MS4A1",
                             "CD79A"),
                split.by = "groups")


## DOT PLOTS
    ## size corresponds to the percentage of cells
    ## expressing the feature in each cluster.
    ## the color represents the average expression level
    DotPlot(pbmc3k.final, features = features) +
      RotatedAxis()
    
    ## SplitDotPlotGG has been replaced with the "split.by" parameter
    DotPlot(pbmc3k.final, 
            features = features,
            split.by = "groups") +
      RotatedAxis()



## HEATMAP
    ## Single cell heatmap of feature expression
    DoHeatmap(subset(pbmc3k.final, downsample = 100),
              features = features,
              size = 3)
    
    ## DoHeatmap shows a grouping bar, splitting the heatmap into groups/clusters
    ## Can be changed with 'group_by' parameter
    DoHeatmap(pbmc3k.final, 
              features = VariableFeatures(pbmc3k.final) [1:100],
              cells = 1:500,
              size = 4,
              angle = 90) + NoLegend()


## DIMP LOTS
    ## DimPlot replaces TSNEPlot, PCAPlot, etc.
    ## in addition, it will plot either "umap", "tsne" or "pca" by default
    ## in that order.
    DimPlot(pbmc3k.final)
    
    pbmc3k.final.no.umap <- pbmc3k.final
    pbmc3k.final.no.umap[["umap"]] <- NULL
    DimPlot(pbmc3k.final.no.umap) + RotatedAxis()


## APPLYING THEMES TO PLOTS
    ## With Seurat, all plotting functions return ggplot2-based plots by default, 
    ## allowing one to easily capture and manipulate plots just like any other 
    ## ggplot2-based plot.
    baseplot <- DimPlot(pbmc3k.final,
                        reduction = "umap")
    ## add custom labels + titles
    baseplot + labs(title = "Clustering of 2,700 PBMCs")
    
    ## there are community-created themes,
    ## they overwrite the default with:
    remotes::install_github('sjessa/ggmin')
    baseplot + ggmin::theme_powerpoint()
    
    ## There are also built-in themes in Seurat
    ## try ?SeuratTheme
    baseplot + DarkTheme()
    
    ## Chain themes together
    baseplot + FontSize(x.title = 20,
                        y.title = 20) +
      NoLegend()
    
    
    
## INTERACTIVE PLOTTING FEATURES

    ## Seurat utilizes R’s plotly graphing library to create interactive plots. 
    ## This interactive plotting feature works with any ggplot2-based scatter plots 
    ## (requires a geom_point layer). To use, simply make a ggplot2-based scatter plot 
    ## (such as DimPlot() or FeaturePlot()) and pass the resulting plot to HoverLocator()
    
        # Include additional data to display alongside cell names by passing in a data frame of
        # information Works well when using FetchData
    plot <- FeaturePlot(pbmc3k.final, features = "MS4A1")
    HoverLocator(plot = plot, 
                     information = FetchData(pbmc3k.final,
                                             vars = c("ident",
                                                      "PC_1",
                                                      "nFeature_RNA")))
    
    ## Another interactive feature provided by Seurat is being able to manually select cells
    ## for further investigation. We have found this particularly useful for small clusters 
    ## that do not always separate using unbiased clustering, but which look tantalizingly distinct. 
    ## You can now select these cells by creating a ggplot2-based scatter plot (such as with DimPlot() 
    ## or FeaturePlot(), and passing the returned plot to CellSelector(). 
    ## CellSelector() will return a vector with the names of the points selected, 
    ## so that you can then set them to a new identity class and perform differential expression.
    
    ## For example, lets pretend that DCs had merged with monocytes in the clustering, 
    ## but we wanted to see what was unique about them based on their position in the tSNE plot.
    
    
    pbmc3k.final <- RenameIdents(pbmc3k.final, 
                                 DC = "CD14+ Mono")
    plot <- DimPlot(pbmc3k.final, 
                    reduction = "umap")
    select.cells <- CellSelector(plot = plot)
    
    ## We can change the identity of these cells to turn them into their own mini-cluster.
    ## select a mini-cluster on the map
    head(select.cells)
    Idents(pbmc3k.final, cells = select.cells) <- "NewCells"
    
    ## Now we can find markers that are specific to the new cells, and find clear DC markers
    newcells.markers <- FindMarkers(pbmc3k.final, 
                                    ident.1 = "NewCells",
                                    ident.2 = "CD14+ Mono",
                                    min.diff.pct = 0.3,
                                    only.pos = T)
    head(newcells.markers)
    
## PLOTTING ACCESSORIES
    ## Along with new functions add interactive functionality to plots, 
    ## Seurat provides new accessory functions for manipulating and combining plots.

    ## LabelClusters and LabelPoints will label the clusters (a coloring variable)
    ## on a ggplot2-based scatter plot
    plot <- DimPlot(pbmc3k.final, 
                    reduction = "pca") + NoLegend()
    LabelClusters(plot = plot, 
                  id = "ident")
    
    ## both functions support 'repel', which will intelligently stagger labels
    ## and draw connection lines from the labels to the points or clusters
    LabelPoints(plot = plot,
                points = TopCells(object = pbmc3k.final[["pca"]]),
                repel = T)
    
    ## Plotting multiple plots was previously achieved with the CombinePlot() function. 
    ## We are deprecating this functionality in favor of the patchwork system. 
    ## Below is a brief demonstration but please see the patchwork package website 
    ## for more details and examples.
    
    plot1 <- DimPlot(pbmc3k.final)
    plot2 <- FeatureScatter(pbmc3k.final, 
                            feature1 = "LYZ",
                            feature2 = "CCL5")
    plot1 + plot2    
    
    ## you can remove the legend from all plots
    (plot1 + plot2) & NoLegend()
    
  