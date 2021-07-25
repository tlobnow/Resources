#install.packages("Seurat")
library(Seurat)
#install.packages("remotes")
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)
library(cowplot)
library(dplyr)
#InstallData("bmcite")
library(bmcite.SeuratData)



#### WNN ANALYSIS OF CITE-SEQ, RNA + ADT #######################################

## CITE-Seq dataset by Stuart, Butler et al, Cell 2019
## it consists of 30.672 scRNA-seq profiles 
## measured alongside a panel of 25 antibodies from Bone Marrow.
## the object contains 2 assays, RNA and ADT (antibody-derived tags).

bm <- LoadData(ds = "bmcite")

## perform pre-processing and dimensional reduction on both assays independently.
## use Std-Normalization or use SCTransform or alternative methods

DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()


DefaultAssay(bm) <- 'ADT'
## use ADT features for dimensional reduction
## we set a dim.reduction name to avoid overwriting
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>%
  RunPCA(reduction.name = 'apca')
## for each cell, we calc. the closest neighbor in the dataset, based on a weighted
## combo of RNA and protein similarities. the cell-specific modality weights
## and multimodal neighbors are calculated in a single fct. (2 min)
## we specify the dimensionality of each modality (~ specifying the # of PCs to 
## include in scRNA-seq clustering), but you can vary these settings to see that 
## small changes have made minimal effect on the overallresults


## Identify multimodal neighbors, will be stored and can be accessed using 'bm'[['weighted.nn']]
## WNN graph can be accessed --> bm[["wknn"]]
## SNN graph used for clustering --> bm[["wsnn"]]
## cell-specific modality weights can be accessed --> bm$RNA.weight
bm <- FindMultiModalNeighbors(bm, 
                              reduction.list = list("pca", "apca"),
                              dims.list = list(1:30, 1:18),
                              modality.weight.name = "RNA.weight")

## We can now use these results for downstream analysis, such as visualization 
## and clustering. For example, we can create a UMAP visualization of the data 
## based on a weighted combination of RNA and protein data. We can also perform 
## graph-based clustering and visualize these results on the UMAP, 
## alongside a set of cell annotations.
bm <- RunUMAP(bm, 
              nn.name = "weighted.nn",
              reduction.name = "wnn.umap",
              reduction.key = "wnnUMAP_")
bm <- FindClusters(bm,
                   graph.name = "wsnn",
                   algorithm = 3,
                   resolution = 2,
                   verbose = F)
p1 <- DimPlot(bm, 
              reduction = 'wnn.umap',
              label = T,
              repel = T,
              label.size = 2.5) + NoLegend()
p2 <- DimPlot(bm, 
              reduction = 'wnn.umap',
              group.by = 'celltype.l2',
              label = T,
              repel = T,
              label.size = 2.5) + NoLegend()
p1 + p2

## We can also compute UMAP visualization based on only the RNA and protein data and compare. 
## We find that the RNA analysis is more informative than the ADT analysis 
## in identifying progenitor states (the ADT panel contains markers for differentiated cells), 
## while the converse is true of T cell states (where the ADT analysis outperforms RNA).
bm <- RunUMAP(bm, 
              reduction = 'pca',
              dims = 1:30,
              assay = 'RNA',
              reduction.name = 'rna.umap',
              reduction.key = 'rna.umap_')
bm <- RunUMAP(bm, 
              reduction = 'apca',
              dims = 1:18,
              assay = 'ADT',
              reduction.name = 'adt.umap',
              reduction.key = 'adt.umap_')

p3 <- DimPlot(bm, 
              reduction = 'rna.umap',
              group.by = 'celltype.l2',
              label = T,
              repel = T,
              label.size = 2.5) + NoLegend()
p4 <- DimPlot(bm, 
              reduction = 'adt.umap',
              group.by = 'celltype.l2',
              label = T,
              repel = T,
              label.size = 2.5) + NoLegend()
p3 + p4


## We can visualize the expression of canonical marker genes and proteins
## on the multimodal UMAP, which can assist in verifying the provided annotations:

p5 <- FeaturePlot(bm,
                  features = c("adt_CD45RA", "adt_CD16", "adt_CD161"),
                  reduction = 'wnn.umap',
                  max.cutoff = 2,
                  cols = c("lightgrey", "darkgreen"),
                  ncol = 3)
p6 <- FeaturePlot(bm,
                  features = c("rna_TRDC", "rna_MPO", "rna_AVP"),
                  reduction = 'wnn.umap',
                  max.cutoff = 3,
                  ncol = 3)
p5 / p6


## Finally we can visualize the modality weights that were learned for each cell.
## Each of the populations with the highest RNA weights represents progenitor cells,
## while the populations with the highest protein weights represent T cells.
## This is in line with our biological expectations, as the antibody-panel 
## does not contain markers that can distinguish between diff. progenitor populations.
VlnPlot(bm, 
        features = "RNA.weight",
        group.by = 'celltype.l2',
        sort = T,
        pt.size = 0.1) + NoLegend()


