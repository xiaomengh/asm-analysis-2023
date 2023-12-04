library(dplyr)
library(Seurat)
library(cowplot)
library(sctransform)
library(RColorBrewer)
library(ggplot2)

sample = "7pretreat3healthy"
#asm = readRDS("7pretreat3healthy_data.rds")

asm.data=Read10X(data.dir = '../../cellranger/AGG7pretreat3healthy/outs/filtered_feature_bc_matrix/')
asm=CreateSeuratObject(counts = asm.data, project="asm", min.cells = 3)
asm[['percent.mt']] = PercentageFeatureSet(object=asm, pattern="^MT-")

asm = subset(x=asm, subset = nFeature_RNA <5000 & percent.mt<15)
asm = SCTransform(asm)

asm = RunPCA(object = asm, features = VariableFeatures(object=asm))
DimPlot(object = asm, reduction="pca")
asm = FindNeighbors(object = asm, dims = 1:30)
asm = FindClusters(object = asm, resolution = 0.5)
asm = RunUMAP(asm, dims = 1:30)
umapplot=paste(sample,"_umap.pdf",sep="")
pdf(umapplot, width=8, height = 6)
DimPlot(object = asm, reduction = "umap", label=TRUE)
dev.off()

savefile = paste(sample,"_data.rds", sep="")
saveRDS(asm, file=savefile)

umapcoordinate = paste(sample, "_umapCoordinate.csv",sep="")
write.csv(asm@reductions$umap@cell.embeddings, umapcoordinate,quote = F)

###### cell markers ######
cellmarkerplot = paste(sample,"_cell_markers.pdf", sep="")
pdf(cellmarkerplot, width = 16, height = 12)
FeaturePlot(object = asm, features = c("MS4A1","CD79A","MZB1", "NKG7","GNLY", "CD3E", "CD8A","CD4",
                                        "CD14","FCGR3A","FCER1A","FCGR3B",
                                        "PPBP","HBB","CD34"),cols=c("grey","red"))
dev.off()

###### plot umap colored by sample #####
idents = read.csv("7pretreat3healthy_idents.csv",header=F)
colnames(idents)=c("barcode","sample","cluster")
asm$sample = idents$sample
asm$cluster = idents$cluster
Idents(asm)= asm$sample
reorder_levels = c("Pt1","Pt2","Pt3","Pt4","CMML1","CMML2","CMML3","H1","H2","H3")
asm@active.ident = factor(x=asm@active.ident,levels=reorder_levels)
pdf(paste(sample,"_umap_groupbysample.pdf",sep=""),width=7,height=6)
DimPlot(asm,reduction ="umap",cols=c("#f8766d","#d89000","#a2a500","#39b600",
                                     "#ff63bc","#e66bf3","#958fff",
                                     "#aed0e8","#71a5d0","#3457a7"))
dev.off()

#######C9 is plasma cells ####
unknown_cells=read.csv("C9C10C11C12.csv",header=T)
c9= unknown_cells[unknown_cells$C9C10C11C12 == "C9",]$Barcode
c10 = unknown_cells[unknown_cells$C9C10C11C12 == "C10",]$Barcode
c11= unknown_cells[unknown_cells$C9C10C11C12 == "C11",]$Barcode 
c12 = unknown_cells[unknown_cells$C9C10C11C12 == "C12",]$Barcode

####### C11
C11 = subset(asm, cells=c11)
pdf("C11_umap.pdf",width=6, height=4)
DimPlot(object=C11, reduction = "umap",pt.size = 2,cols=c("#f8766d","#d89000","#a2a500","#39b600",
                                                          "#ff63bc","#e66bf3","#958fff",
                                                          "#aed0e8","#71a5d0","#3457a7"))
dev.off()

pdf("C11_releventmarkers.pdf",width=20,height=12)
plist = FeaturePlot(object = C11, features = c("CCR3", "CLC","HDC","MS4A2",
                                               "KIT","TPSAB1","CPA3","FCER1A",
                                               "LAIR1","ITGA4","IL3RA","CD34"),cols=c("grey","red"),combine=F,pt.size = 2)
plist = lapply(plist, function(g){
  g +xlim(c(1,2.5))+ylim(c(-4.2,-2))
})
CombinePlots(plist,ncol=4)
dev.off()


####### C10
C10 = subset(asm, cells=c10)
pdf("C10_umap.pdf",width=6, height=4)
DimPlot(object=C10, reduction = "umap",pt.size=2,cols=c("#f8766d","#d89000","#39b600",
                                                        "#ff63bc",
                                                        "#aed0e8","#3457a7"))
dev.off()
pdf("C10_releventmarkers.pdf",width=20,height=12)
plist = FeaturePlot(object = C10, features = c("CCR3", "CLC","HDC","MS4A2",
                                               "KIT","TPSAB1","CPA3","FCER1A",
                                               "LAIR1","ITGA4","IL5RA","IL3RA"),cols=c("grey","red"),combine=F,pt.size = 2)
plist = lapply(plist, function(g){
  g +xlim(c(-0.95,-0.7))+ylim(c(-4.7,-4.45))
})
CombinePlots(plist,ncol=4)
dev.off()

####### C12
C12 = subset(asm, cells=c12)
pdf("C12_umap.pdf",width=6, height=4)
DimPlot(object=C12, reduction = "umap",pt.size=2,cols=c("#f8766d","#d89000",
                                                        "#958fff"))
dev.off()

pdf("C12_releventmarkers.pdf",width=10,height=12)
plist = FeaturePlot(object = C12, features = c("TUBB2B","COL1A1","COL1A2","COL3A1","TIMP1","SPP1"),cols=c("grey","red"),combine=F,pt.size=2)
plist = lapply(plist, function(g){
  g +xlim(c(2,2.8))+ylim(c(-2.3,-1.5))
})
CombinePlots(plist,ncol=2)
dev.off()

###### Neutrophil ####
Idents(asm) = asm$cluster
neu=subset(asm,idents=c("4","17"))

pdf("neutrophil_bySample.pdf",width=30,height=4)
DimPlot(neu,split.by="sample",pt.size=2)+ylim(-13,-7)
dev.off()
pdf("neutrophil_granulocyte_markers.pdf",width=11,height=5)
plist = FeaturePlot(object=neu,features=c("ITGAM","FCGR3B","CEACAM3","CEACAM1","CEACAM8"),cols=c("grey","red"),pt.size = 0.5,combine=FALSE)
plist = lapply(plist, function(g){
  g +ylim(c(-13,-7))
})
CombinePlots(plist,ncol=3)
dev.off()
