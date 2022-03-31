library(dplyr)
library(Seurat)
library(patchwork)
getwd()
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

## STANDARD PRE-PROCESSING WORKFLOW
    ## The steps below encompass the standard pre-processing workflow for scRNA-seq data 
    ## in Seurat. These represent the selection and filtration of cells based on QC metrics, 
    ## data normalization and scaling, and the detection of highly variable features.


## QC AND SELECTING CELLS FOR FURTHER ANALYSIS
    ## Seurat allows you to easily explore QC metrics and filter cells based on 
    ## any user-defined criteria. A few QC metrics commonly used by the community include:

          ## The number of unique genes detected in each cell.
              ## Low-quality cells or empty droplets will often have very few genes
              ## Cell doublets or multiplets may exhibit an aberrantly high gene count
          ## Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
          ## The percentage of reads that map to the mitochondrial genome
              ## Low-quality / dying cells often exhibit extensive mitochondrial contamination
              ## We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
              ## We use the set of all genes starting with MT- as a set of mitochondrial genes


    ## The [[ operator can add columns to object metadata. 
    ## This is a great place to stash QC Stats!
        pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc,
                                                 pattern = "^MT-")
    
        ## Where are QC matrices stored in Seurat??
        ## in the example, we visualize QC metrics, and use these to filter cells
            ## we filter cells with unique feature counts over 2500 or under 200
            ## we filter cells that have >5% mitochondrial counts
        ## we will do this using a Violin Plot
        VlnPlot(pbmc, features = c("nFeature_RNA",
                                   "nCount_RNA",
                                   "percent.mt"),
                ncol = 3)
    
    ## FeatureScatter is typically used to visualize feature-feature relationships,
    ## but can be used for anything calculated by the object
        ## i.e. columns in object metadata, PC scores, etc.
        
    plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    plot1 + plot2
        
    pbmc <- subset(pbmc, 
                   subset = 
                     nFeature_RNA > 200 & 
                     nFeature_RNA < 2500 &
                     percent.mt < 5)
    
    
## NORMALIZING THE DATA      
    ## After removing unwanted cells from the dataset, the next step is to 
    ## normalize the data. By default, we employ a global-scaling normalization 
    ## method “LogNormalize” that normalizes the feature expression measurements 
    ## for each cell by the total expression, multiplies this by a scale factor 
    ## (10,000 by default), and log-transforms the result. 
    ## Normalized values are stored in pbmc[["RNA"]]@data.
    pbmc <- NormalizeData(pbmc,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)
    
    ## For clarity, in the previous line of code (and in future commands), 
    ## we provide the default values for certain parameters in the function call.
    ## However, this isn't required and the same behavior can also be achieved with:
    pbmc <- NormalizeData(pbmc)
    
## IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION)
    ## We next calculate a subset of features that exhibit 
    ## high cell-to-cell variation in the dataset 
    ## (i.e, they are highly expressed in some cells, and lowly expressed in others). 
    ## We and others have found that focusing on these genes in downstream analysis 
    ## helps to highlight biological signal in single-cell datasets.
    
    ## Our procedure in Seurat is described in detail here, 
    ## and improves on previous versions by directly modeling the mean-variance 
    ## relationship inherent in single-cell data, and is implemented in the 
    ## FindVariableFeatures() function. By default, we return 2,000 features per dataset. 
    ## These will be used in downstream analysis, like PCA.

    pbmc <- FindVariableFeatures(pbmc,
                                 selection.method = "vst",
                                 nfeatures = 2000)
    ## identify the 10 most highly variable genes:
    top10 <- head(VariableFeatures(pbmc), 10)
    
    ## plot variable features with and without labels
    plot1 <- VariableFeaturePlot(pbmc)
    plot2 <- LabelPoints(plot = plot1,
                         points = top10,
                         repel = T)
    plot1 + plot2
        

## SCALING THE DATA
    ## Now we will apply a linear transformation ("scaling") that is a standard
    ## pre-processing step prior to dimensional reduction techniques like PCA.
    ## The ScaleData() function:
        ## Shifts the expression of each gene, so that the mean expression across cells is 0
        ## Scales the expression of each gene, so that the variance across cells is 1:
              ## This step gives equal weight in downstream analyses, 
              ## so that highly-expressed genes do not dominate
        ## The results are stored in pbmc[["RNA"]]@scale.data
    
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc,
                      features = all.genes)

    ## This step can take pretty long..
        ## Scaling is an essential step in the Seurat workflow, but only on genes 
        ## that will be used as input to PCA. Therefore, the default in ScaleData() 
        ## is only to perform scaling on the previously identified variable features 
        ## (2,000 by default). To do this, omit the features argument in the previous 
        ## function call, i.e.: pbmc <- ScaleData(pbmc)
        ## Your PCA and clustering results will be unaffected. However, 
        ## Seurat heatmaps (produced as shown below with DoHeatmap()) require genes 
        ## in the heatmap to be scaled, to make sure highly-expressed genes don’t 
        ## dominate the heatmap. To make sure we don’t leave any genes out of the heatmap later, 
        ## we are scaling all genes in this tutorial.
    
    ## How to remove unwanted sources of variation, as in Seurat v2:
        ## In Seurat v2 we also use the ScaleData() function to remove unwanted 
        ## sources of variation from a single-cell dataset. 
        ## For example, we could ‘regress out’ heterogeneity associated with (for example) 
        ## cell cycle stage, or mitochondrial contamination. 
        ## These features are still supported in ScaleData() in Seurat v3, i.e.:
          ## pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
        ## However, particularly for advanced users who would like to use this functionality, 
        ## we strongly recommend the use of our new normalization workflow, SCTransform(). 
        ## The method is described in our paper, with a separate vignette using Seurat v3 here. 
        ## As with ScaleData(), the function SCTransform() also includes a vars.to.regress parameter.
    
## PERFORM LINEAR DIMENSIONAL REDUCTION
    ## next we will perform PCA on the scaled data.
    ## by default, only the previously determined variable features are used as input,
    ## but you can define this using the 'features' argument, 
    ## if you wish to choose a different subset.
    pbmc <- RunPCA(pbmc,
                   features = VariableFeatures(object = pbmc))
    
    ## Seurat provides several useful ways of visualizing both cells and features
    ## that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
    ## Let's look at different ways to visualize PCA results!
    print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
    
    VizDimLoadings(pbmc, 
                   dims = 1:2,
                   reduction = "pca")
    DimPlot(pbmc, 
            reduction = "pca")
    
    ## In particular DimHeatmap() allows for easy exploration of the primary sources
    ## of heterogeneity in a dataset, and can be useful when trying to decide
    ## which PCs to include for further downstream analyses.
    ## Both cells and features are ordered according to their PCA scores.
    ## Setting cells to a number plots 'extreme' cells on both ends of the spectrum,
    ## which dramatically speeds plotting for large datasets.
    ## Though clearly a supervised analysis, we find this to be a valuable tool
    ## for exploring correlated feature sets.
    
    DimHeatmap(pbmc,
               dims = 1,
               cells = 500,
               balanced = T)

    DimHeatmap(pbmc,
               dims = 1:15,
               cells = 500,
               balanced = T)
    
    
## DETERMINE THE 'DIMENSIONALITY' OF THE DATASET
    ## To overcome the extensive technical noise in any single feature for scRNA-seq data, 
    ## Seurat clusters cells based on their PCA scores, with each PC essentially 
    ## representing a ‘metafeature’ that combines information across a correlated feature set. 
    ## The top principal components therefore represent a robust compression of the dataset. 
    ## However, how many components should we choose to include? 10? 20? 100?
    
    ## In Macosko et al, we implemented a resampling test inspired by the JackStraw procedure. 
    ## We randomly permute a subset of the data (1% by default) and rerun PCA, 
    ## constructing a ‘null distribution’ of feature scores, and repeat this procedure. 
    ## We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.

      ## NOTE: This process can take a long time for big datasets, comment out for expediency. More
      ## approximate techniques such as those implemented in ElbowPlot() can be used to reduce
      ## computation time
    pbmc <- JackStraw(pbmc, num.replicate = 100)
    pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
      
      ## The JackStrawPlot()  function provides a visualization tool for comparing 
      ## the distribution of p-values for each PC with a uniform distribution (dashed line). 
      ## ‘Significant’ PCs will show a strong enrichment of features with low p-values 
      ## (solid curve above the dashed line). In this case it appears that there is a 
      ## sharp drop-off in significance after the first 10-12 PCs.
    JackStrawPlot(pbmc, dims = 1:15)

    ## An alternative heuristic method generates an 'Elbow Plot'
    ## It's a ranking of principle components based on the percentage of variance
    ## explained by each one. In this example, we can observe an 'elbow' around PC9/10,
    ## suggesting that the majority of true signal is captured in the first 10 PCs
    ElbowPlot(pbmc)

    ## Identifying the true dimensionality of a dataset – can be challenging/uncertain for the user. 
    ## We therefore suggest these three approaches to consider. 
    ## The first is more supervised, exploring PCs to determine relevant sources 
    ## of heterogeneity, and could be used in conjunction with GSEA for example.
    ## The second implements a statistical test based on a random null model, 
    ## but is time-consuming for large datasets, and may not return a clear PC cutoff. 
    ## The third is a heuristic that is commonly used, and can be calculated instantly. 
    ## In this example, all three approaches yielded similar results, 
    ## but we might have been justified in choosing anything between PC 7-12 as a cutoff.
    
    ## We chose 10 here, but encourage to consider the following points:
      ## Dendritic cell and NK aficionados may recognize that genes strongly 
      ## associated with PCs 12 and 13 define rare immune subsets 
      ## (i.e. MZB1 is a marker for plasmacytoid DCs). 
      ## However, these groups are so rare, they are difficult to distinguish from 
      ## background noise for a dataset of this size without prior knowledge.
    
      ## We encourage users to repeat downstream analyses with a different number 
      ## of PCs (10, 15, or even 50!). As you will observe, the results often 
      ## do not differ dramatically.
    
      ## We advise users to err on the higher side when choosing this parameter. 
      ## For example, performing downstream analyses with only 5 PCs 
      ## does significantly and adversely affect results.
      
## CLUSTERING THE CELLS
    ## Seurat v3 applies a graph-based clustering approach, building upon initial 
    ## strategies in (Macosko et al). Importantly, the distance metric which drives 
    ## the clustering analysis (based on previously identified PCs) remains the same. 
    ## However, our approach to partitioning the cellular distance matrix into 
    ## clusters has dramatically improved. Our approach was heavily inspired by 
    ## recent manuscripts which applied graph-based clustering approaches to 
    ## scRNA-seq data [SNN-Cliq, Xu and Su, Bioinformatics, 2015] and 
    ## CyTOF data [PhenoGraph, Levine et al., Cell, 2015]. Briefly, 
    ## these methods embed cells in a graph structure - for example a K-nearest 
    ## neighbor (KNN) graph, with edges drawn between cells with similar feature 
    ## expression patterns, and then attempt to partition this graph into highly 
    ## interconnected ‘quasi-cliques’ or ‘communities’.
    
    ## As in PhenoGraph, we first construct a KNN graph based on the euclidean 
    ## distance in PCA space, and refine the edge weights between any two cells 
    ## based on the shared overlap in their local neighborhoods (Jaccard similarity). 
    ## This step is performed using the FindNeighbors() function, and takes as 
    ## input the previously defined dimensionality of the dataset (first 10 PCs).
    
    ## To cluster the cells, we next apply modularity optimization techniques 
    ## such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., 
    ## Journal of Statistical Mechanics], to iteratively group cells together, 
    ## with the goal of optimizing the standard modularity function. 
    ## The FindClusters() function implements this procedure, and contains a 
    ## resolution parameter that sets the ‘granularity’ of the downstream clustering, 
    ## with increased values leading to a greater number of clusters. 
    ## We find that setting this parameter between 0.4-1.2 typically returns 
    ## good results for single-cell datasets of around 3K cells. 
    ## Optimal resolution often increases for larger datasets. 
    ## The clusters can be found using the Idents() function.
    
    pbmc <- FindNeighbors(pbmc, dims = 1:10)
    pbmc <- FindClusters(pbmc, resolution = 0.5)
    
    ## let's look at the Cluster IDs of the first 5 cells:
    head(Idents(pbmc), 5)

## RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/TSNE)
    ## Seurat offers several non-linear dimensional reduction techniques, 
    ## such as tSNE and UMAP, to visualize and explore these datasets. 
    ## The goal of these algorithms is to learn the underlying manifold of the 
    ## data in order to place similar cells together in low-dimensional space. 
    ## Cells within the graph-based clusters determined above should co-localize 
    ## on these dimension reduction plots. As input to the UMAP and tSNE, 
    ## we suggest using the same PCs as input to the clustering analysis.
    
    pbmc <- RunUMAP(pbmc, dims = 1:10)
    ## you can set label = T or use the labelClusters() fct. 
    ## to help label individual clusters
    DimPlot(pbmc, reduction = "umap")
    

## FINDING DIFFERENTIALLY EXPRESSED FEATURES (CLUSTER BIOMARKERS)
    ## Seurat can help you find markers that define clusters via differential expression. 
    ## By default, it identifies positive and negative markers of a single cluster 
    ## (specified in ident.1), compared to all other cells. FindAllMarkers() 
    ## automates this process for all clusters, but you can also test groups of 
    ## clusters vs. each other, or against all cells.
    
    ## The min.pct argument requires a feature to be detected at a minimum percentage 
    ## in either of the two groups of cells, and the thresh.test argument requires 
    ## a feature to be differentially expressed (on average) by some amount between 
    ## the two groups. You can set both of these to 0, but with a dramatic increase 
    ## in time - since this will test a large number of features that are unlikely 
    ## to be highly discriminatory. As another option to speed up these computations, 
    ## max.cells.per.ident can be set. This will downsample each identity class 
    ## to have no more cells than whatever this is set to. While there is generally 
    ## going to be a loss in power, the speed increases can be significant and 
    ## the most highly differentially expressed features will likely still rise 
    ## to the top.
    
    ## find all markers of cluster 2
    cluster2.markers <- FindMarkers(pbmc,
                                    ident.1 = 2,
                                    min.pct = 0.25)
    head(cluster2.markers)

    ## find all markers distinguishing cluster 5 from clusters 0 and 3
    cluster5.markers <- FindMarkers(pbmc,
                                    ident.1 = 5,
                                    ident.2 = c(0,3),
                                    min.pct = 0.25)
    head(cluster5.markers, 
         n = 5)

    ## find markers for every cluster compared to all remaining cells,
    ## report only the positive ones:
    pbmc.markers <- FindAllMarkers(pbmc,
                                   only.pos = T,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25)
    pbmc.markers %>%
      group_by(cluster) %>%
      top_n(n = 2,
            wt = avg_log2FC)
    ## Seurat has several tests for differential expression which can be set 
    ## with the test.use parameter (see our DE vignette for details). 
    ## For example, the ROC test returns the ‘classification power’ 
    ## for any individual marker (ranging from 0 - random, to 1 - perfect).
    cluster0.markers <- FindMarkers(pbmc, 
                                    ident.1 = 0,
                                    logfc.threshold = 0.25,
                                    test.use = "roc",
                                    only.pos = T)
    
    ## We include several tools for visualizing marker expression. VlnPlot() 
    ## (shows expression probability distributions across clusters), 
    ## and FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) 
    ## are our most commonly used visualizations. 
    ## We also suggest exploring RidgePlot(), CellScatter(), and DotPlot() 
    ## as additional methods to view your dataset.
    
    VlnPlot(pbmc,
            features = c("MS4A1",
                         "CD79A"))

    ## You can plot raw counts too:
    VlnPlot(pbmc,
            features = c("NKG7",
                         "PF4"),
            slot = "counts",
            log = T)

    FeaturePlot(pbmc,
                features = c("MS4A1", 
                             "GNLY", 
                             "CD3E", 
                             "CD14", 
                             "FCER1A", 
                             "FCGR3A", 
                             "LYZ", 
                             "PPBP",
                             "CD8A"))    
    
    ## DoHeatmap() generates an expression heatmap for given cells and features.
    ## In this case, we are plotting the top20 markers (or all markes if <20) for each cluster
    pbmc.markers %>%
      group_by(cluster) %>%
      top_n(n = 10,
            wt = avg_log2FC) -> top10
    DoHeatmap(pbmc,
              features = top10$gene) + NoLegend()
    
## ASSIGNING CELL TYPE IDENTITY TO CLUSTERS
    ## fortunately in this dataset, we can use canonical markers to easily
    ## match the unbiased clustering to known cell types!
    # Cluster ID	Markers	       Cell Type
    #   0	    IL7R, CCR7	      Naive CD4+ T
    #   1	    CD14, LYZ	        CD14+ Mono
    #   2	    IL7R, S100A4	    Memory CD4+
    #   3	    MS4A1	            B
    #   4	    CD8A	            CD8+ T
    #   5	    FCGR3A, MS4A7	    FCGR3A+ Mono
    #   6	    GNLY, NKG7	      NK
    #   7	    FCER1A, CST3	    DC
    #   8	    PPBP	            Platelet
    
    new.cluster.ids <- c("Naive CD4 T", 
                         "CD14+ Mono", 
                         "Memory CD4 T", 
                         "B", 
                         "CD8 T", 
                         "FCGR3A+ Mono",
                         "NK", 
                         "DC", 
                         "Platelet")
    names(new.cluster.ids) <- levels(pbmc)
    pbmc <- RenameIdents(pbmc,
                         new.cluster.ids)
    DimPlot(pbmc,
            reduction = "umap",
            label = T,
            pt.size = 0.5) + NoLegend()
    
    ## you can save your data as:
    #saveRDS(pbmc, file = "pbmc3k_final.rds")
    
    
    
    
    
    