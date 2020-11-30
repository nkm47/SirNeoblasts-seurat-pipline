##############################
library(dplyr)
library(Seurat)

###  import data and remove mito RNA , rRNA    
##############################################################
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


### Normalizing the data, Identification of highly variable features, and Scaling      
##############################################################

# log Normalize
sirneoblast <- NormalizeData(sirneoblast, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
sirneoblast <- FindVariableFeatures(sirneoblast, selection.method = "vst", nfeatures = 2000)

# Scaling the data  
all.genes <- rownames(sirneoblast)
sirneoblast <- ScaleData(sirneoblast, features = all.genes)


###  Perform linear dimensional reduction  and Determine the ‘dimensionality’ of the dataset
##############################################################

sirneoblast <- RunPCA(sirneoblast, features = VariableFeatures(object = sirneoblast),npcs = 50)
print(sirneoblast[["pca"]], dims = 1:30, nfeatures = 5)

# JackStraw
sirneoblast <- JackStraw(sirneoblast, num.replicate = 100,dims = 40)
print(JackStrawPlot(sirneoblast, dims = 1:30))

# Elbow plot方法
sirneoblast <- ScoreJackStraw(sirneoblast, dims = 1:30)
print(ElbowPlot(sirneoblast,ndims=30))


### Cluster the cells       
##############################################################

pca_num<- 10  
sirneoblast <- FindNeighbors(sirneoblast, dims = 1:pca_num)    
res_cluster<- 0.6  
sirneoblast <- FindClusters(sirneoblast, resolution = res_cluster)
head(Idents(sirneoblast), 5)

### Run non-linear dimensional reduction (tSNE)  
##############################################################
set.seed(7)
sirneoblast <- RunTSNE(sirneoblast, dims = 1:pca_num)

# tSNE plot
tsne_name_png<- paste("8.cluster-",res_cluster,"label-tsne-pca-",pca_num,".png",sep="")
png(tsne_name_png,width = 1100,height = 1000, res=72*3)
print(DimPlot(sirneoblast, reduction = "tsne",label.size = 3, label = TRUE,pt.size = 0.08))
dev.off()
rm(tsne_name_png)


###  Finding differentially expressed features  
##############################################################
sirneoblast.markers <- FindAllMarkers(sirneoblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
maker_name <- paste("SirNeoblasts-pca_", pca_num, "_resolution_", res_cluster, "all_markers.csv", sep = "")
write.csv(sirneoblast.markers, maker_name)
