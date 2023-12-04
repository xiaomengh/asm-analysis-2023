library(dplyr)
library(Seurat)
library(cowplot)
library(sctransform)
library(ggplot2)

sample = "4pts_3H"
#asm = readRDS("4pts_3H_data.rds")

asm.data=Read10X(data.dir = '../../cellranger/4pts_3H/outs/count/filtered_feature_bc_matrix/')
asm=CreateSeuratObject(counts = asm.data, project="asm", min.cells = 3)
asm[['percent.mt']] = PercentageFeatureSet(object=asm, pattern="^MT-")
FeatureScatter(object = asm, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
FeatureScatter(object = asm, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(object = asm, feature1 = "nCount_RNA", feature2 = "percent.mt")

asm = subset(x=asm, subset = nFeature_RNA <5000 & percent.mt<15)
asm = SCTransform(asm)

asm = RunPCA(object = asm, features = VariableFeatures(object=asm))
DimPlot(object = asm, reduction="pca")

####### Cluster the cells #############
asm = FindNeighbors(object = asm, dims = 1:30)
asm = FindClusters(object = asm, resolution = 0.5)
asm = RunUMAP(asm, dims = 1:30)
umapplot=paste(sample,"_umap.pdf",sep="")
pdf(umapplot, width=8, height = 6)
DimPlot(object = asm, reduction = "umap", label=TRUE)
dev.off()

###### save umap coordinate #####
umapcoordinate = paste(sample, "_umapCoordinate.csv",sep="")
write.csv(asm@reductions$umap@cell.embeddings, umapcoordinate,quote = F)

###### save data #####
savefile = paste(sample,"_data.rds", sep="")
saveRDS(asm, file=savefile)

###### cell markers ######
cellmarkerplot = paste(sample,"_cell_markers.pdf", sep="")
pdf(cellmarkerplot, width = 16, height = 12)
FeaturePlot(object = asm, features = c("MS4A1","CD79A","MZB1", "NKG7","GNLY", "CD3E", "CD8A","CD4",
                                        "CD14","FCGR3A","FCER1A","FCGR3B",
                                        "PPBP","HBB","CD34"),cols=c("grey","red"))
dev.off()

###### plot umap colored by pt and sample #####
idents = read.csv("4pts_3H_idents_new.csv",header=F)
colnames(idents)=c("barcode","sample","tp","cluster")
asm$sample = idents$sample
asm$tp = idents$tp
asm$cluster=idents$cluster
asm$sample_tp = paste(idents$sample,idents$tp,sep="_")

Idents(asm) = asm$sample_tp
reorder_levels=c("Pt1_T1","Pt1_T2","Pt1_T3",
                 "Pt2_T1","Pt2_T2","Pt2_T3",
                 "Pt3_T1","Pt3_T2","Pt3_T3",
                 "Pt4_T1","Pt4_T2","Pt4_T3",
                 "H1_T0","H2_T0","H3_T0")
asm@active.ident=factor(x=asm@active.ident,levels=reorder_levels)
pdf(paste(sample,"_umap_groupbySample.pdf",sep=""), width = 5, height = 4)
DimPlot(object=asm,reduction="umap")
dev.off()

####### cleanup monocytes data, remove CD16+ mono and DC #####
Idents(asm) = cluster
mono = subset(asm,idents=c("2","12","20"))
DimPlot(mono,reduction="umap",split.by="sample")
mono=RunPCA(obj=mono,features = VariableFeatures(object=mono))
mono=RunUMAP(mono,dim=1:10)
mono = FindNeighbors(object = mono, dims = 1:10)
mono = FindClusters(object = mono, resolution = 0.5)
DimPlot(object=mono, reduction = "umap", label=TRUE)
FeaturePlot(object = mono, features = c("CD14","FCER1A","FCGR3A", "FCGR3B"))
mono_p = subset(mono,idents=c("0","1","2","3","4","5"))
mono_p=RunPCA(obj=mono_p,features = VariableFeatures(object=mono))
mono_p=RunUMAP(mono_p,dim=1:20)
mono_p = FindNeighbors(object = mono_p, dims = 1:20)
mono_p = FindClusters(object = mono_p, resolution = 0.3)
DimPlot(object=mono_p, reduction = "umap", label=TRUE)
mono_viable= subset(mono_p,idents=c("0","1","2"))
saveRDS(mono_viable,"4pts_3H_viable_mono_data.rds")

### longitudinal cell umap
cellUmapPlotbySample=function(cell_data,pt_id,colors=c("#ff6c67","#bda007","#0dbf26","#aed0e8","#71a5d0","#3457a7")){
  sample_names = rownames(table(Idents(cell_data)))
  ident_need = c(paste(pt_id,"_T1",sep=""),
                 paste(pt_id,"_T2",sep=""),
                 paste(pt_id,"_T3",sep=""),
                 "H1_T0","H2_T0","H3_T0")
  ident_sel= ident_need[ident_need %in% sample_names]
  pt_data = subset(cell_data,idents=ident_sel)
  print(table(Idents(pt_data)))
  p1 = DimPlot(pt_data,reduction="umap",cols=colors,pt.size=2)
  return(p1)
}
# monocytes
mono_viable=readRDS("4pts_3H_viable_mono_data.rds")
Idents(mono_viable) = mono_viable$sample_tp
new.cluster.ids=c("Pt2_T1","Pt2_T2","Pt2_T3",
                  "Pt3_T1","Pt3_T2","Pt3_T3",
                  "Pt4_T1","Pt4_T2","Pt4_T3",
                  "Pt1_T1","Pt1_T2","Pt1_T3",
                  "H1_T0","H2_T0","H3_T0"
)
names(new.cluster.ids)=levels(mono_viable)
mono_viable = RenameIdents(mono_viable,new.cluster.ids)
reorder_levels=c("Pt1_T1","Pt1_T2","Pt1_T3",
                 "Pt2_T1","Pt2_T2","Pt2_T3",
                 "Pt3_T1","Pt3_T2","Pt3_T3",
                 "Pt4_T1","Pt4_T2","Pt4_T3",
                 "H1_T0","H2_T0","H3_T0")
mono_viable@active.ident=factor(x=mono_viable@active.ident,levels=reorder_levels)
write.csv(Idents(mono_viable),"mono_viable.idents.csv",quote=F)
pdf("Pt1_3H_mono_umap.pdf",width=5,height=4)
cellUmapPlotbySample(mono_viable,"Pt1")
dev.off()
pdf("Pt2_3H_mono_umap.pdf",width=5,height=4)
cellUmapPlotbySample(mono_viable,"Pt2")
dev.off()
pdf("Pt3_3H_mono_umap.pdf",width=5,height=4)
cellUmapPlotbySample(mono_viable,"Pt3")
dev.off()
pdf("Pt4_3H_mono_umap.pdf",width=5,height=4)
cellUmapPlotbySample(mono_viable,"Pt4")
dev.off()

#neutrophil
Idents(asm) = asm$cluster
neu = subset(asm,idents = c("3","4"))
Idents(neu)=neu$sample_tp
table(Idents(neu))
pdf("Pt1_3H_neutrophil_umap.pdf",width=5,height=4)
cellUmapPlotbySample(neu, "Pt1",colors=c("#ff6c67","#bda007","#0dbf26","#aed0e8","#3457a7"))+xlim(c(-5,8))+ylim(c(6,12))
dev.off()
pdf("Pt2_3H_neutrophil_umap.pdf",width=5,height=4)
cellUmapPlotbySample(neu, "Pt2",colors=c("#ff6c67","#bda007","#0dbf26","#aed0e8","#3457a7"))+xlim(c(-5,8))+ylim(c(6,12))
dev.off()
pdf("Pt3_3H_neutrophil_umap.pdf",width=5,height=4)
cellUmapPlotbySample(neu, "Pt3",colors=c("#ff6c67","#0dbf26","#aed0e8","#3457a7"))+xlim(c(-5,8))+ylim(c(6,12))
dev.off()
pdf("Pt4_3H_neutrophil_umap.pdf",width=5,height=4)
cellUmapPlotbySample(neu, "Pt4",colors=c("#ff6c67","#bda007","#0dbf26","#aed0e8","#3457a7"))+xlim(c(-5,8))+ylim(c(6,12))
dev.off()

ima_neu = subset(asm,idents = c("16"))
Idents(ima_neu)=ima_neu$sample_tp
table(Idents(ima_neu))
pdf("Pt1_3H_immature_neutrophil_umap.pdf",width=5,height=4)
cellUmapPlotbySample(ima_neu, "Pt1",colors=c("#ff6c67","#bda007","#0dbf26"))
dev.off()
pdf("Pt2_3H_immature_neutrophil_umap.pdf",width=5,height=4)
cellUmapPlotbySample(ima_neu, "Pt2",colors=c("#ff6c67","#bda007","#0dbf26"))
dev.off()
pdf("Pt3_3H_immature_neutrophil_umap.pdf",width=5,height=4)
cellUmapPlotbySample(ima_neu, "Pt3",colors=c("#ff6c67","#0dbf26"))
dev.off()
pdf("Pt4_3H_immature_neutrophil_umap.pdf",width=5,height=4)
cellUmapPlotbySample(ima_neu, "Pt4",colors=c("#ff6c67","#bda007"))
dev.off()

allneu = subset(asm,idents = c("3","4","16"))
Idents(allneu)=allneu$sample_tp
table(Idents(allneu))
write.csv(Idents(allneu),"allneu.idents.csv",quote=F)

#######
unknow = read.csv("Ba_C11_manual.csv",header=T)
ba=unknow[unknow$Ba_C11=="Ba",]$Barcode
c11 = unknow[unknow$Ba_C11=="C11",]$Barcode

Ba = subset(asm, cells=ba)
Idents(Ba) = Ba$sample_tp
table(Idents(Ba))
pdf("Pt1_3H_Ba_umap.pdf",width=5,height=4)
cellUmapPlotbySample(Ba, "Pt1",colors=c("#ff6c67","#bda007","#0dbf26","#aed0e8","#3457a7"))
dev.off()
pdf("Pt2_3H_Ba_umap.pdf",width=5,height=4)
cellUmapPlotbySample(Ba, "Pt2",colors=c("#ff6c67","#aed0e8","#3457a7"))
dev.off()
pdf("Pt3_3H_Ba_umap.pdf",width=5,height=4)
cellUmapPlotbySample(Ba, "Pt3",colors=c("#0dbf26","#aed0e8","#3457a7"))
dev.off()
pdf("Pt4_3H_Ba_umap.pdf",width=5,height=4)
cellUmapPlotbySample(Ba, "Pt4",colors=c("#ff6c67","#bda007","#0dbf26","#aed0e8","#3457a7"))
dev.off()

C11 = subset(asm, cells=c11)
Idents(C11) = C11$sample_tp
table(Idents(C11))
pdf("Pt1_3H_C11_umap.pdf",width=5,height=4)
cellUmapPlotbySample(C11, "Pt1",colors=c("#ff6c67","#bda007","#0dbf26"))+xlim(c(-5,-2))+ylim(c(2,3.5))
dev.off()
pdf("Pt2_3H_C11_umap.pdf",width=5,height=4)
cellUmapPlotbySample(C11, "Pt2",colors=c("#ff6c67","#bda007","#0dbf26"))+xlim(c(-5,-2))+ylim(c(2,3.5))
dev.off()
pdf("Pt3_3H_C11_umap.pdf",width=5,height=4)
cellUmapPlotbySample(C11, "Pt3",colors=c("#ff6c67","#bda007","#0dbf26"))+xlim(c(-5,-2))+ylim(c(2,3.5))
dev.off()
pdf("Pt4_3H_C11_umap.pdf",width=5,height=4)
cellUmapPlotbySample(C11, "Pt4",colors=c("#ff6c67","#bda007","#0dbf26"))+xlim(c(-5,-2))+ylim(c(2,3.5))
dev.off()

####### NK, T and B cells
# NK
Idents(asm)=asm$seurat_clusters
lym=subset(asm,idents = 10)
Idents(lym)=lym$sample_tp
table(Idents(lym))
lym_celltype="NK"
pdf(paste("Pt1_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt1")+xlim(c(-1,6))+ylim(c(-13,-8))
dev.off()
pdf(paste("Pt2_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt2")+xlim(c(-1,6))+ylim(c(-13,-8))
dev.off()
pdf(paste("Pt3_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt3")+xlim(c(-1,6))+ylim(c(-13,-8))
dev.off()
pdf(paste("Pt4_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt4",colors=c("#ff6c67","#0dbf26","#aed0e8","#71a5d0","#3457a7"))+xlim(c(-1,6))+ylim(c(-13,-8))
dev.off()

# T
nkt=read.csv("T_NKcells_manual.csv",header=T)
tcells = nkt$Barcode[!nkt$Barcode %in% idents$barcode[idents$cluster==10]]
lym=subset(asm,cells=tcells)
Idents(lym)=lym$sample_tp
table(Idents(lym))
lym_celltype="T"
pdf(paste("Pt1_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt1")+xlim(c(-15,6))+ylim(c(-15,0))
dev.off()
pdf(paste("Pt2_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt2")+xlim(c(-15,6))+ylim(c(-15,0))
dev.off()
pdf(paste("Pt3_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt3")+xlim(c(-15,6))+ylim(c(-15,0))
dev.off()
pdf(paste("Pt4_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt4")+xlim(c(-15,6))+ylim(c(-15,0))
dev.off()

# B
Idents(asm)=asm$seurat_clusters
lym=subset(asm,idents = 14)
Idents(lym)=lym$sample_tp
table(Idents(lym))
lym_celltype="B"
# plot
pdf(paste("Pt1_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt1")+xlim(c(-15,-10))+ylim(c(5,8))
dev.off()
pdf(paste("Pt2_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt2")+xlim(c(-15,-10))+ylim(c(5,8))
dev.off()
pdf(paste("Pt3_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt3")+xlim(c(-15,-10))+ylim(c(5,8))
dev.off()
pdf(paste("Pt4_3H_",lym_celltype,"_umap.pdf",sep=""),width=5,height=4)
cellUmapPlotbySample(lym, "Pt4")+xlim(c(-15,-10))+ylim(c(5,8))
dev.off()


########## calculate cell percentage #####
celllineage = read.csv("cell_lineage.csv",header=T)
Idents(asm) = asm$sample_tp
reorder_levels=c("Pt1_T1","Pt1_T2","Pt1_T3",
                 "Pt2_T1","Pt2_T2","Pt2_T3",
                 "Pt3_T1","Pt3_T2","Pt3_T3",
                 "Pt4_T1","Pt4_T2","Pt4_T3",
                 "H1_T0","H2_T0","H3_T0")
asm@active.ident=factor(x=asm@active.ident,levels=reorder_levels)

nucleated_cell = subset(asm,cell=celllineage[celllineage$cell.lineage!="plt_RBC",]$Barcode)
table(Idents(nucleated_cell))

myeloid=subset(asm, cell=celllineage[celllineage$cell.lineage=="myeloid_cells",]$Barcode)
table(Idents(myeloid))
#########################
