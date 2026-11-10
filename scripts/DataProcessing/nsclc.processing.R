########### NSCLC processing


#nsclc preprocessing

nsclc<-NormalizeData(nsclc)
nsclc <- FindVariableFeatures(nsclc, selection.method = "vst", nfeatures = 2000)
nsclc<-ScaleData(nsclc)
nsclc<-RunPCA(nsclc)
nsclc <- RunHarmony(nsclc, 
                        group.by.vars = c( "orig.ident", "day"),
                        reduction = "pca", 
                        assay.use = "SCT", 
                        reduction.save = "harmony")
nsclc <- RunUMAP(nsclc, reduction = "harmony", dims = 1:20)
nsclc <- FindNeighbors(nsclc, reduction = "harmony", dims = 1:20)
nsclc <- FindClusters(nsclc.obj, resolution = 1.64)

nsclc <- FindAllMarkers(nsclc, min.pct = 0.25, logfc.threshold = 0.25, only.pos=T)

## annotated based on DEGs
Idents(nsclc) <- "All.Cell.Subtypes"
nsclc <- RenameIdents(nsclc,
                      "CD8 Teff(1)" = "CD8 T",
                      "CD8 Teff(2)" = "CD8 T",
                      "CD8 Teff(3)" = "CD8 T",
                      "CD8 Teff(4)" = "CD8 T",
                      "CD8 Teff(5)" = "CD8 T",
                      "CD8 Trm" = "CD8 T",
                      "CD8 Texh" = "CD8 T",
                      "CD4-CD8 IFN" = "CD8 T",
                      # "CD4-CD8 IFN" = "CD8 T IFN",
                      "CD4 naive-like" = "CD4 Th",
                      "CD4 Tcm" = "CD4 Th",
                      "CD4 Treg" = "CD4 Treg",
                      "CD4 Tfh-like" = "CD4 Tfh",
                      "T apoptotic" = "T undefined",
                      "T unspecified" = "T undefined",
                      "Plasma cells" = "Plasma B cell",
                      "B cells" = "B cell",
                      "CD16+ NK" = "CD16+ NK",
                      "XCL1 NK" = "NK",
                      "ILC" = "ILC",
                      "Mast cell" = "Mast cell",
                      "CD14+ Mono" = "CD14+ Monocyte",
                      "CD16+ Mono" = "CD16+ Monocyte",
                      "Antimicrobial MoDCs" = "MoDC",
                      "SLC11A1 MoDCs" = "MoDC",
                      "Ag-presenting cDC2" = "cDC2",
                      "IFN cDC2" = "cDC2",
                      "cDC1" = "cDC1",
                      "migDCs" = "migDC",
                      "Macrophages" = "Macrophage APOE",
                      "Myeloid unspecified(1)" = "Macrophage LYZ",
                      "Myeloid unspecified(2)" = "Macrophage LYZ",
                      "pDC(1)" = "pDC",
                      "pDC(2)" = "pDC",
                      "pDC(3)" = "pDC",
                      "stromal" = "Stromal cell",
                      "Proliferating(1)" = "Proliferating",
                      "Proliferating(2)" = "Proliferating", # in figure proliferating 1
                      "Cycling" = "Cycling" # in figure 'proliferating 2
                      
)
nsclc@meta.data$All.Subsets <- nsclc@active.ident