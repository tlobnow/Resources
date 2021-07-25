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
    
    ## Interactive Plotting Features







