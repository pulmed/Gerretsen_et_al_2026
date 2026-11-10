### split T cells from nsclc object

Idents(nsclc) <- "ann_1"
T.tcr <- subset(nsclc.tcr, idents = c("CD8","CD4 Tem", "CD4 Treg", "CD4 Tfh"))

## add NK and cycling cells with TCRs
test1 <- subset(nsclc.tcr, idents = c("cycling","NK"))
test1 <- subset(x = test1, subset = (Frequency >= "1"))
T.tcr <- merge(T.tcr, test1, add.cell.ids = c("T", "o"))

DefaultAssay(T.tcr) <- 'RNA'
T.tcr<-FindVariableFeatures(T.tcr)
T.tcr<-NormalizeData(T.tcr)
T.tcr<-ScaleData(T.tcr)
T.tcr<-RunPCA(T.tcr)
T.tcr <- RunHarmony(T.tcr, 
                    group.by.vars = c( "orig.ident", "day"),
                    reduction = "pca", 
                    assay.use = "SCT", 
                    reduction.save = "harmony")
T.tcr<- RunUMAP(T.tcr, reduction = "harmony", assay = "SCT", dims = 1:20)
T.tcr <- FindNeighbors(object = T.tcr, dims = 1:25, reduction = "harmony")
T.tcr <- FindClusters(T.tcr, resolution = c(0.3, 1.0)) 
DimPlot(T.tcr, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T, pt.size = 1)


# Check markers
FeaturePlot(T.tcr, c('CD3D', 'CD3E'))
FeaturePlot(T.tcr, c('CD4', 'CD8A', 'CD8B','FOXP3'))
FeaturePlot(T.tcr, c('SELL', 'CD69', 'IL7R','TCF7'))
VlnPlot(T.tcr, c( "CD8A", "CD8B"))
VlnPlot(T.tcr, c( "adt_CD4-A0045", "adt_CD4-A0072"), pt.size=0)
FeaturePlot(T.tcr, c( "adt_CD4-A0045", "adt_CD4-A0072"))

T.AM <-FindAllMarkers(T.cell) # performed on 0.3

Idents(T.tcr) <- "RNA_snn_res.1"
T.tcr <- RenameIdents(T.tcr,
                      #T4
                      "11" = "C1",
                      "0" = "C1",
                      "1" = "C2",
                      "5"= "C3",
                      "14" = "C4",
                      "3" = "C5",
                      "17" = "C5",
                      "8" = "C6",
                      "15" = "C6",
                      "6" = "C7",
                      #T8
                      "19" = "C8",
                      "7" = "C9",
                      "9" = "C10",
                      "10" = "C11",
                      "12" = "C12",
                      "4" = "C13",
                      "2" = "C14",
                      "13" = "C15",
                      "16" = "C16",
                      #T4-T8
                      "18" = "C17",
                      "20" = "C18",
                      "21" = "C19")

T.tcr@meta.data$annot_T5 <- Idents(T.tcr)
Idents(T.tcr) <- "annot_T5"
T4.tcr <- subset(T.tcr, idents = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C17", "C18", "C19"))
T8.tcr <- subset(T.tcr, idents = c("C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16"))

## split into T8 and T4
#Idents(T.tcr) <- "ann_T1"
#T4 <- subset(T.tcr, idents = c("CD4.TEMRA" , "CD4.N", "CD4.CM", "CD4.EM", "CD4.Tfh/EM", "CD4.Treg"))
#T8 <- subset(T.tcr, idents = c("CD8.TEMRA", "CD8.EM" ,  "CD8.RM", "CD8.cycling", "CD8.proliferating"))

## Identify CD8 cells in T4 and CD4 cells in T8 object
T4n <- subset(T4, idents = c('0', '7'), invert = T)
T8v <- subset(T4, idents = c('0', '7'))
T8n <- subset(T8, idents = c('0','1','9'), invert = T)
T4v <- subset(T8, idents = c('0','1','9'))

# Pull cell names
cells1 <- colnames(T4n)
cells2 <- colnames(T4v)

# Identify shared cells between both objects
shared_cells <- intersect(x = cells1, y = cells2)

# Create Obj1 minus shared cells
T4v <- subset(T4v, cells = shared_cells, invert = TRUE)

# Pull cell names
cells1 <- colnames(T8n)
cells2 <- colnames(T8v)

# Identify shared cells between both objects
shared_cells <- intersect(x = cells1, y = cells2)

# Create Obj1 minus shared cells
T8v <- subset(T8v, cells = shared_cells, invert = TRUE)

# Merge objects
T4 <- merge(T4n, T4v)
T8 <- merge(T8n, T8v)


saveRDS(T8, "T8.rds")
saveRDS(T4, "T4.rds")



###################################################################################
## T8 processing
DefaultAssay(T8) <- 'RNA'
T8<-FindVariableFeatures(T8)
T8<-NormalizeData(T8)
T8<-ScaleData(T8)
T8<-RunPCA(T8)
T8 <- RunHarmony(T8, 
                 group.by.vars = c( "orig.ident", "day"),
                 reduction = "pca", 
                 assay.use = "SCT", 
                 reduction.save = "harmony")
ElbowPlot(T8, ndims = 30, reduction = "harmony")
T8<- RunUMAP(T8, reduction = "harmony", assay = "SCT", dims = 1:27) # 27 dims based on elbowplot
ElbowPlot(T8, ndims = 30, reduction = "harmony")
T8 <- FindNeighbors(object = T8, dims = 1:18, reduction = "harmony") # changed dims based on elbowplot
T8 <- FindClusters(T8, resolution = 1.3)
genes.use <- grep(pattern = "^RPS|^RPL|^HSP|^MT-|^MALAT1|^DNAJ|^LINC0", ##
                  rownames(T8),
                  value=TRUE, invert=TRUE) 
Idents(T8) <- "RNA_snn_res.1.3"
T_AM8 <- FindAllMarkers(T8, min.pct = 0.25, features = genes.use) %>% Add_Pct_Diff()
Idents(T8) <- "RNA_snn_res.1.3"
T8<- RenameIdents(T8,
                  "9"= "CD8-GZMH", # GNLY GZMH - cytotoxic / transitional
                  "19" = "CD8-LTB", #early naive?
                  "0" = "CD8-LTB",  # early naive?
                  "5" = "CD8-SELL", # CCR7 NOSIP LEF1 TCF7
                  "2" = "CD8-KLRC1", # ook ZNF683, IL7R, HOPX, 
                  "13" = "CD8-KLF2", # CD8 Tn CCR7 LEF1 TCF7 SELL NOSIP
                  "3" = "CD8-TNFSF9", #41BB, GZMK
                  "8" = "CD8-TNFSF9", #41BB, NR4A1, EGR1
                  "7" = "CD8-TNFAIP3", # samen met 10
                  "10" = "CD8-TNFAIP3", # samen met 7
                  "1" = "CD8-GZMK", # 
                  
                  "6"= "CD8-ZNF683", # CXCL13 (hobit)
                  "4"= "CD8-RGS1", # CTLA4 CXCL13 CCL5
                  "11"= "CD8-KLRC2", # XCL1, XCL2, ZNF683 HOPX
                  "14"= "CD8-KLRC2",# XCL1, XCL2, ZNF683 HOPX
                  
                  "16"= "CD8-CTLA4", # TNFSF4
                  "12"= "CD8-CXCL13", #GZMK CXCR4
                  "15"= "CD8-MKI67", #EOMES-CD27
                  "18"= "CD8-FABP5", #-B1,A,B5,DQ, GZMK
                  "17" = "CD8-KLRB1"
)
saveRDS(T8, "T8.rds")


########################################################
## T4 processing
DefaultAssay(T4) <- 'RNA'
T4<-FindVariableFeatures(T4)
T4<-NormalizeData(T4)
T4<-ScaleData(T4)
T4<-RunPCA(T4)
T4 <- RunHarmony(T4, 
                 group.by.vars = c( "orig.ident", "day"),
                 reduction = "pca", 
                 assay.use = "SCT", 
                 reduction.save = "harmony")
T4<- RunUMAP(T4, reduction = "harmony", assay = "SCT", dims = 1:20)
T4 <- FindNeighbors(object = T4, dims = 1:25, reduction = "harmony")
T4 <- FindClusters(T4, resolution = 0.5)
Idents(T4) <- "RNA_snn_res.0.5"
DimPlot(T4, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T, pt.size = 1) + DimPlot(T4, reduction = "umap", group.by = "ann_2", label = T, pt.size = 1) 
genes.use <- grep(pattern = "^RPS|^RPL|^HSP|^MT-|^MALAT1|^MKI67|^STMN1", ##
                  rownames(T4),
                  value=TRUE, invert=TRUE) #get list of non-ribosomal genes 
T4_AM <- FindAllMarkers(T4, min.pct = 0.25, features = genes.use) %>% Add_Pct_Diff()

Idents(T4) <- "cluster_annotation"

T4_cd4$cluster_annotation <- Idents(T4_cd4)

T4 <- RenameIdents(T4,
                   "Tearly act"="Tearly act", 
                   "Tcm1"="Tcentral mem(1)",
                   "Tcm2"="Tcentral mem(2)",
                   "T HSP"="T HSP",
                   "Tcytotox1"="Tcytotox(1)",
                   "Tcytotox2"="Tcytotox(2)",
                   "Tcytotox3"="Tcytotox(3)",
                   "Tem mix" = "CD4CD8",
                   "Treg1"="Treg(1)",
                   "Treg1"="Treg(2)",
                   "Texh-like"="Texh-like",
                   "Tfh"="Tfh")

saveRDS(T4, "T4.rds")



