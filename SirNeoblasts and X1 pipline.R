##############################
library(dplyr)
library(Seurat)

###  import data and remove mito RNA , rRNA    
##############################################################
#X1
x1.data<- read.csv("GSE107873_X1_cells_10x_exprs.csv",header = T,row.names = 1)
x1 <- CreateSeuratObject(counts = x1.data, project = "X1")
dim(x1)
x1[["percent.mt"]] <- PercentageFeatureSet(x1, features = c("SMED30031308", "SMED30034623" ,"SMED30019959", 
                                   "SMED30010375" ,"SMED30012212" ,"SMED30008414" ,
                                   "SMED30026531","SMED30001799" ,"SMED30005014" ,
                                   "SMED30013300","SMED30017284", "SMED30024106", 
                                   "SMED30012974","SMED30022128","SMED30003672" ,
                                   "SMED30024888","SMED30024665","SMED30000702",
                                   "SMED30032663" ,"SMED30027845","SMED30032887"))
x1 <- subset(x1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 70)
dim(x1)

#print(FeatureScatter(x1, feature1 = "percent.mt", feature2 = "nFeature_RNA"))
count_data<-as.data.frame(x1@assays$RNA@counts)
count_data<-count_data[-match(c("SMED30031308", "SMED30034623" ,"SMED30019959", 
                                   "SMED30010375" ,"SMED30012212" ,"SMED30008414" ,
                                   "SMED30026531","SMED30001799" ,"SMED30005014" ,
                                   "SMED30013300","SMED30017284", "SMED30024106", 
                                   "SMED30012974","SMED30022128","SMED30003672" ,
                                   "SMED30024888","SMED30024665","SMED30000702",
                                   "SMED30032663" ,"SMED30027845","SMED30032887"),rownames(count_data)),]
x1<- CreateSeuratObject(counts = count_data, project = "x1-an")
dim(x1)


#SirNeoblast
sirneoblast<- Read10X(data.dir = "./sir_filtered_gene_bc_matrices")
sirneoblast<- CreateSeuratObject(counts = sirneoblast, project = "sirNeoblast", min.cells = 3, min.features = 200)
dim(sirneoblast)
sirneoblast[["percent.mt"]] <- PercentageFeatureSet(sirneoblast, features = c("SMED30031308", "SMED30034623" ,"SMED30019959", 
                                   "SMED30010375" ,"SMED30012212" ,"SMED30008414" ,
                                   "SMED30026531","SMED30001799" ,"SMED30005014" ,
                                   "SMED30013300","SMED30017284", "SMED30024106", 
                                   "SMED30012974","SMED30022128","SMED30003672" ,
                                   "SMED30024888","SMED30024665","SMED30000702",
                                   "SMED30032663" ,"SMED30027845","SMED30032887"))
sirneoblast <- subset(sirneoblast, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 70)

count_data<-as.data.frame(sirneoblast@assays$RNA@counts)
count_data<-count_data[-match(c("SMED30031308", "SMED30034623" ,"SMED30019959", 
                                   "SMED30010375" ,"SMED30012212" ,"SMED30008414" ,
                                   "SMED30026531","SMED30001799" ,"SMED30005014" ,
                                   "SMED30013300","SMED30017284", "SMED30024106", 
                                   "SMED30012974","SMED30022128","SMED30003672" ,
                                   "SMED30024888","SMED30024665","SMED30000702",
                                   "SMED30032663" ,"SMED30027845","SMED30032887"),rownames(count_data)),]


# data combined
##############################################################
all_sample<- list(sirneoblast,x1)

for (i in 1:length(all_sample)) {
  all_sample[[i]] <- NormalizeData(all_sample[[i]], verbose = FALSE)
  all_sample[[i]] <- FindVariableFeatures(all_sample[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

sdata<- FindIntegrationAnchors(object.list = all_sample, dims = 1:10)
data.combined <- IntegrateData(anchorset = sdata, dims = 1:10)
#Perform an integrated analysis
DefaultAssay(data.combined ) <- "integrated"


# Run the standard workflow for clustering
##############################################################
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:10)
#resolution
data.combined <- FindClusters(data.combined, resolution = 0.6 )

### Run non-linear dimensional reduction (tSNE)  
##############################################################
#tsne
set.seed(7)
data.combined <- RunTSNE(data.combined, reduction = "pca", dims = 1:10)
p1 <- DimPlot(data.combined, reduction = "tsne", group.by = "orig.ident",pt.size = 0.5)  
p2 <- DimPlot(data.combined, reduction = "tsne", label = F,pt.size = 0.08)
png("tsne-two-combined.png",width = 3100,height = 1600, res=72*3)
print(plot_grid(p1, p2))
dev.off()
png("tsne-groupby-orig-ident.png",width = 1100,height = 800, res=72*2)
print(p1)
dev.off()
png("3-tsne-plot.png",width = 1100,height = 1000, res=72*3)
print(p2)
dev.off()
p3<-DimPlot(data.combined, reduction = "tsne", split.by = "orig.ident",pt.size = 0.5)
p3
png("tsne-split.png",width = 1600,height = 800, res=72*2)
print(p3)
dev.off()
p4<-DimPlot(data.combined, reduction = "tsne", split.by = "orig.ident",label = TRUE)
png("tsne-split-label.png",width = 2100,height = 1600, res=72*2)
print(p4)
dev.off()


# Finding differentially expressed features
##############################################################
data.combined.markers <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
maker_name <- paste("SirNeoblasts.and.X1-pca_", pca_num, "_resolution_", res_cluster, "all_markers.csv", sep = "")
write.csv(data.combined.markers, maker_name)