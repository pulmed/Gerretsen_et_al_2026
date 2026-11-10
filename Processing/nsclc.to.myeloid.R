## Myeloid subclustering

# to get better and more representative separation of myeloid cells
nsclc$Myeloid <- ifelse(Idents(nsclc) %in% c("macrophages", "class.mono", "non.class.mono", "cDC2", "cDC1", "migDC") & UMAP[, "UMAP_1"] > 8 & UMAP[, "UMAP_1"] < 15 & UMAP[, "UMAP_2"] > -3 & UMAP[, "UMAP_2"] < 5, "Yes", "No")
nsclc.Myeloid <- subset(nsclc, Myeloid == "Yes")
nsclc.Myeloid$ann_3 <- factor(nsclc.Myeloid$ann_3, levels = levels(nsclc.Myeloid$ann_3)[table(nsclc.Myeloid$ann_3) > 0])
# somehow the nsclc.Myeloid cannot be used directly on my end, thus create a new SeuratObject
Myeloid <- nsclc.Myeloid[["RNA"]]@counts
Myeloid.meta.data <- nsclc.Myeloid@meta.data[, c("orig.ident", "ann_3")]
Myeloid <- CreateSeuratObject(counts = Myeloid, meta.data = Myeloid.meta.data)
DefaultAssay(Myeloid) <- "RNA"
Myeloid[["RNA"]] <- split(Myeloid[["RNA"]], f = Myeloid$orig.ident)
Myeloid <- NormalizeData(Myeloid)
Myeloid <- FindVariableFeatures(Myeloid)
Myeloid <- ScaleData(Myeloid)
Myeloid <- RunPCA(Myeloid)
Myeloid <- FindNeighbors(Myeloid, dims = 1:30, reduction = "pca")
Myeloid <- FindClusters(Myeloid, resolution = 2, cluster.name = "unintegrated_clusters")
Myeloid <- RunUMAP(Myeloid, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(Myeloid, reduction = "umap.unintegrated", group.by = c("orig.ident", "ann_3"))
# cluster mainly by orig.ident, thus, do integration
# Seurat 5 integration, to remove batch effect
Myeloid <- IntegrateLayers(
  object = Myeloid, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
Myeloid <- FindNeighbors(Myeloid, reduction = "integrated.cca", dims = 1:50)
Myeloid <- FindClusters(Myeloid, resolution = 2, cluster.name = "cca_clusters")
Myeloid <- RunUMAP(Myeloid, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")
DimPlot(Myeloid, reduction = "umap.cca", group.by = c("orig.ident", "ann_3", "cca_clusters"), combine = TRUE, label.size = 2)
# orig.ident integrated, myeloid cells more separated, just looking at the DimPlot with existing ann_3

# resolution too high, starting from the resolution that result in similar number of cluster with ann_3, which is cca_clusters.0.3
Myeloid.join <- JoinLayers(Myeloid)
Myeloid.join <- FindClusters(Myeloid.join, cluster.name = "cca_clusters.0.3", resolution = 0.3)
DimPlot(Myeloid.join, group.by = "cca_clusters.0.3", reduction = "umap.cca") | DimPlot(Myeloid.join, group.by = "orig.ident", reduction = "umap.cca")
# well integrated

table(Myeloid.join$orig.ident, Myeloid.join$ann_3)
# some samples are more enriched for some clusters, like nsclcE009_reseq and nsclcE014_reseq for non.class.mono
table(Myeloid.join$orig.ident, Myeloid.join$cca_clusters.0.3)
# cca_clusters.0.3 similar to ann_3
VlnPlot(Myeloid.join, features = c("nCount_RNA", "nFeature_RNA"), group.by = "orig.ident")
# nsclcE001_reseq more lymphocytes, thus lower nCount_RNA and lower nFeature_RNA
cca_clusters.0.3.DEGs <- FindAllMarkers(Myeloid.join, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, recorrect_umi = FALSE) # recorrect_umi = FALSE since it is subset
cca_clusters.0.3.DEGs$gene[cca_clusters.0.3.DEGs$cluster == 0]# Exhausted lymphocytes, aligns with mmc6, Lymphocytes exhausted/cell-cycle, from Sade_Feldman, table s5 from https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(24)00035-1?uuid=uuid%3A94e8fbcb-b8ae-4fc0-b609-5bc20c1d6a06
cca_clusters.0.3.DEGs$gene[cca_clusters.0.3.DEGs$cluster == 1]# DC2, this and following naming are based on classic markers and also mmc8 (somehow no reference paper, but the markers are reliable)
cca_clusters.0.3.DEGs$gene[cca_clusters.0.3.DEGs$cluster == 2]# CD14+ Mono (class.mono)
cca_clusters.0.3.DEGs$gene[cca_clusters.0.3.DEGs$cluster == 3]# Mac
cca_clusters.0.3.DEGs$gene[cca_clusters.0.3.DEGs$cluster == 4]# cycling
FeaturePlot(Myeloid.join, features = c("C5AR1", "CD3E"), reduction = "umap.cca") # cca_clusters.0.3 == 4 is mainly lymphocytes that express CD3E
FeaturePlot(Myeloid.join, features = c("STMN1", "MKI67"), reduction = "umap.cca") # STMN1 is more specific than MKI67, which should make it "cycling STMN1", I will leave it to you and call it Cycling instead here
cca_clusters.0.3.DEGs$gene[cca_clusters.0.3.DEGs$cluster == 5]# DC1
cca_clusters.0.3.DEGs$gene[cca_clusters.0.3.DEGs$cluster == 6]# mig.DCs
cca_clusters.0.3.DEGs$gene[cca_clusters.0.3.DEGs$cluster == 7]# CD16+ Mono (non.class.mono)
cca_clusters.0.3.DEGs$gene[cca_clusters.0.3.DEGs$cluster == 8]# Doublets, mix of cells, IG heavy, RP heavy
# saveRDS(cca_clusters.0.3.DEGs, "/Volumes/FindAllMarkers/DEGs/cca_clusters.0.3.norecorrect.rds")

Myeloid.join$Cell.Types <- factor(mapvalues(Myeloid.join$Cell.Subtypes, from = c("Macrophages", "CD14+ Mono", "CD16+ Mono", "cDC1", "Ag-presenting cDC2", "IFN cDC2", "Antimicrobial MoDCs", "SLC11A1 MoDCs", "migDCs", "Lym-Mye Doublets", "Exhausted Lymphocytes", "Cycling"), to = c("Macrophages", "CD14+ Mono", "CD16+ Mono", "cDC1", "cDC2", "cDC2", "MoDCs", "MoDCs", "migDCs", "Lym-Mye Doublets", "Exhausted Lymphocytes", "Cycling")), 
                                  levels = c("Macrophages", "CD14+ Mono", "CD16+ Mono", "cDC1", "cDC2", "MoDCs", "migDCs", "Lym-Mye Doublets", "Exhausted Lymphocytes", "Cycling"))


# saveRDS(Myeloid.join, "/Volumes/SeuratObject/Myeloid.join.rds")
# saveRDS(Myeloid.join@meta.data, "/Volumes/meta.data/Myeloid.join_meta.data.rds")