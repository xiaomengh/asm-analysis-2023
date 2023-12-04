plotCellVariantMatrix = function(data,mainlabel,cells,variant_range=NA){
  if (!is.na(variant_range)) {
    matrix = data[variant_range,][,colnames(data) %in% cells]
    matrix = cbind(data[variant_range,]$variant,matrix)
  } else {
    matrix = data[,colnames(data) %in% cells]
    matrix = cbind(data$variant,matrix)
  }
  #print(head(matrix[,c(1:5)]))
  colnames(matrix) = c("variant", 1:(ncol(matrix)-1))
  matrix$variant = c(1:nrow(matrix))
  #print(head(matrix[,c(1:5)]))
  df_heatmap=reshape2::melt(matrix,id.vars="variant")
  #print(head(df_heatmap))
  df_heatmap$variant = as.numeric(df_heatmap$variant)
  df_heatmap$variable = as.numeric(df_heatmap$variable)
  df_heatmap$value=as.numeric(df_heatmap$value)
  #print(head(df_heatmap))
  b = df_heatmap[df_heatmap$value==0,]
  c = df_heatmap[df_heatmap$value>0,]
  plot(b$variable, b$variant,pch=16, cex=1,col="#F5DEB3",ylim=c(1,nrow(matrix)),xlim=c(1,ncol(matrix)),
       xlab="Cell index",ylab="Variant index",main=mainlabel)
  points(c$variable, c$variant,pch=16, cex=1,col="#80008080",ylim=c(1,nrow(matrix)),xlim=c(1,ncol(matrix)))
}
plotCellVariantBySubcloneAndByCellType = function(cells, celltype,is.all){
  print(celltype)
  sc1_cell_index = colnames(data) %in% cells & (colnames(data) %in% sc1)
  sc2_cell_index = colnames(data) %in% cells & (colnames(data) %in% sc2)
  sc3_cell_index = colnames(data) %in% cells & (colnames(data) %in% sc3)
  normal_index = colnames(data) %in% cells & (colnames(data) %in% normal)
  unassign_index = colnames(data) %in% cells & (colnames(data) %in% unassign)
  sc1_data = data[,sc1_cell_index]
  sc2_data = data[,sc2_cell_index]
  sc3_data = data[,sc3_cell_index]
  print("sc1 cell number")
  print(sum(sc1_cell_index))
  print("sc2 cell number")
  print(sum(sc2_cell_index))
  print("sc3 cell number")
  print(sum(sc3_cell_index))
  print("normal cells")
  print(sum(normal_index))
  print("unassigned cells")
  print(sum(unassign_index))
  normal_data = data[,normal_index]
  unassign_data = data[,unassign_index]
  
  if (is.all){
    new_data = cbind(data$variant,sc1_data,sc2_data,sc3_data,normal_data,unassign_data)
  } else {
    new_data = cbind(data$variant,sc1_data,sc2_data,sc3_data,normal_data)
  }
  colnames(new_data)[1]="variants"
  plotCellVariantMatrix(new_data,mainlabel=paste("genotype for ",celltype,sep=""),colnames(new_data)[2:ncol(new_data)])
  abline(h = 989.5)
  abline(h = 1829.5)
  #abline(h = 2473.5)
  #abline(h = 2524.5)
  abline(v = sum(sc1_cell_index)+0.5)
  abline(v = sum(sc1_cell_index)+sum(sc2_cell_index)+0.5)
  abline(v = sum(sc1_cell_index)+sum(sc2_cell_index)+sum(sc3_cell_index)+0.5) #when include normal
  if (is.all) {
    abline(v = sum(sc1_cell_index)+sum(sc2_cell_index)+sum(sc3_cell_index)+0.5)
    abline(v = sum(sc1_cell_index)+sum(sc2_cell_index)+sum(sc3_cell_index)+ncol(normal_data)+0.5)
  }
}
plotCellVariantBySubcloneSortByQualAndByCellType = function(cells, celltype,is.all,assignment){
  order = order(assignment$QUAL,decreasing=T)
  ordered_data = data[,c(1,order+1)]
  sc1_cell_index = colnames(ordered_data) %in% cells & (colnames(ordered_data) %in% sc1)
  sc2_cell_index = colnames(ordered_data) %in% cells & (colnames(ordered_data) %in% sc2)
  sc3_cell_index = colnames(ordered_data) %in% cells & (colnames(ordered_data) %in% sc3)
  normal_index = colnames(ordered_data) %in% cells & (colnames(ordered_data) %in% normal)
  unassign_index = colnames(ordered_data) %in% cells & (colnames(ordered_data) %in% unassign)
  sc1_data = ordered_data[,sc1_cell_index]
  sc2_data = ordered_data[,sc2_cell_index]
  sc3_data = ordered_data[,sc3_cell_index]
  print("sc1 cell number")
  print(sum(sc1_cell_index))
  print("sc2 cell number")
  print(sum(sc2_cell_index))
  print("sc3 cell number")
  print(sum(sc3_cell_index))
  print("normal cells")
  print(sum(normal_index))
  print("unassigned cells")
  print(sum(unassign_index))
  normal_data = ordered_data[,normal_index]
  unassign_data = ordered_data[,unassign_index]
  
  if (is.all){
    new_data = cbind(ordered_data$variant,sc1_data,sc2_data,sc3_data,normal_data,unassign_data)
  } else {
    new_data = cbind(ordered_data$variant,sc1_data,sc2_data,sc3_data,normal_data)
  }
  colnames(new_data)[1]="variants"
  plotCellVariantMatrix(new_data,mainlabel=paste("genotype for ",celltype,sep=""),colnames(new_data)[2:ncol(new_data)])
  abline(h = 989.5)
  abline(h = 1829.5)
  #abline(h = 2473.5)
  #abline(h = 2524.5)
  abline(v = sum(sc1_cell_index)+0.5)
  abline(v = sum(sc1_cell_index)+sum(sc2_cell_index)+0.5)
  abline(v = sum(sc1_cell_index)+sum(sc2_cell_index)+sum(sc3_cell_index)+0.5) #when include normal
  if (is.all) {
    abline(v = sum(sc1_cell_index)+sum(sc2_cell_index)+sum(sc3_cell_index)+0.5)
    abline(v = sum(sc1_cell_index)+sum(sc2_cell_index)+sum(sc3_cell_index)+ncol(normal_data)+0.5)
  }
}
fillInQual=function(assignment){
  c =matrix(0,nrow=nrow(assignment),ncol=ncol(assignment)-4)
  for(i in 1:(ncol(assignment)-4)){
    a = strsplit(assignment[,4+i],split=":")
    b = matrix(unlist(a),ncol=5,byrow=T)
    c[,i]=as.double(b[,5])
  }
  assignment$QUAL = apply(c,1,max)
  return(assignment)
}
plotUmap = function(xlim,ylim){
  plot(umap$X,umap$Y,pch=16, cex=1,col="white",xlim=xlim, ylim=ylim,xlab="UMAP_1",ylab="UMAP_2")
  points(umap[rownames(umap) %in% normal,][,1],umap[rownames(umap) %in% normal,][,2], pch=16, cex=1,col="#808080")
  points(umap[rownames(umap) %in% sc1,][,1],umap[rownames(umap) %in% sc1,][,2], pch=16, cex=1,col="#7570B4")
  points(umap[rownames(umap) %in% sc2,][,1],umap[rownames(umap) %in% sc2,][,2], pch=16, cex=1,col="#D95F01")
  points(umap[rownames(umap) %in% sc3,][,1],umap[rownames(umap) %in% sc3,][,2], pch=16, cex=1,col="#1C9E76")
}


idents_all=read.csv("../../../4pts_3tp_scRNAseq_analysis/cell_lineage.csv",header=T)
idents_all$Barcode = gsub("-",".",idents_all$Barcode)


########## plot myeloid cells Figure 5B #########
pdf("myeloid_cells_subclone.pdf",width=13,height=4)
par(mfrow=c(1,3))
for (sampleid in c("T1","T2","T3")){
  data=read.table(paste(sampleid,".combined.scaf.tsv",sep=""),header=T)
  #only look at sc1,2,3 variant
  data=data[c(1:2473),]
  if(sampleid=="T1"){
    colnames(data) = gsub(".1",".7",colnames(data))
    assign=read.table(paste(sampleid,".pre.assign",sep=""),sep='\t',header=T)
    assign = fillInQual(assign)
    assign[assign$QUAL<7,]$ASIG = "UNASSIGN"
    assign$Barcode = gsub("-1",".7",assign$Barcode)
  }
  if(sampleid=="T2"){
    colnames(data) = gsub(".1",".8",colnames(data))
    assign=read.table(paste(sampleid,".post.assign",sep=""),sep='\t',header=T)
    assign = fillInQual(assign)
    assign[assign$QUAL<8,]$ASIG = "UNASSIGN"
    assign$Barcode = gsub("-1",".8",assign$Barcode)
  }
  if(sampleid=="T3"){
    colnames(data) = gsub(".1",".9",colnames(data))
    assign=read.table(paste(sampleid,".post.assign",sep=""),sep='\t',header=T)
    assign = fillInQual(assign)
    assign[assign$QUAL<8,]$ASIG = "UNASSIGN"
    assign$Barcode = gsub("-1",".9",assign$Barcode)
  }
  normal = assign[assign$ASIG=='normal',]$Barcode
  sc1 = assign[assign$ASIG=='SC1',]$Barcode
  sc2 = assign[assign$ASIG=='SC2',]$Barcode
  sc3 = assign[assign$ASIG=='SC3',]$Barcode
  unassign = assign[assign$ASIG=='UNASSIGN',]$Barcode
  is.all=F
  cells = idents_all[idents_all$cell.lineage=="myeloid_cells",]$Barcode
  plotCellVariantBySubcloneAndByCellType(cells=cells,celltype="myeloid cells",is.all=is.all)
}
dev.off()

########## plot myeloid cells Suppl. Fig 7 #########
pdf("myeloid_with_unassign.pdf",width=13,height=4)
par(mfrow=c(1,3))
for (sampleid in c("T1","T2","T3")){
  data=read.table(paste(sampleid,".combined.scaf.tsv",sep=""),header=T)
  #only look at sc1,2,3 variant
  data=data[c(1:2473),]
  if(sampleid=="T1"){
    colnames(data) = gsub(".1",".7",colnames(data))
    assign=read.table(paste(sampleid,".pre.assign",sep=""),sep='\t',header=T)
    assign = fillInQual(assign)
    assign[assign$QUAL<7,]$ASIG = "UNASSIGN"
    assign$Barcode = gsub("-1",".7",assign$Barcode)
  }
  if(sampleid=="T2"){
    colnames(data) = gsub(".1",".8",colnames(data))
    assign=read.table(paste(sampleid,".post.assign",sep=""),sep='\t',header=T)
    assign = fillInQual(assign)
    assign[assign$QUAL<8,]$ASIG = "UNASSIGN"
    assign$Barcode = gsub("-1",".8",assign$Barcode)
  }
  if(sampleid=="T3"){
    colnames(data) = gsub(".1",".9",colnames(data))
    assign=read.table(paste(sampleid,".post.assign",sep=""),sep='\t',header=T)
    assign = fillInQual(assign)
    assign[assign$QUAL<8,]$ASIG = "UNASSIGN"
    assign$Barcode = gsub("-1",".9",assign$Barcode)
  }
  normal = assign[assign$ASIG=='normal',]$Barcode
  sc1 = assign[assign$ASIG=='SC1',]$Barcode
  sc2 = assign[assign$ASIG=='SC2',]$Barcode
  sc3 = assign[assign$ASIG=='SC3',]$Barcode
  unassign = assign[assign$ASIG=='UNASSIGN',]$Barcode
  is.all=T
  cells = idents_all[idents_all$cell.lineage=="myeloid_cells",]$Barcode
  plotCellVariantBySubcloneSortByQualAndByCellType(cells=cells,celltype="myeloid cells",is.all=is.all,assign)
}
dev.off()

########## plot lymphocytes from three timepoint together
pdf("lymphocyte_with_unassign.pdf",width=13,height=4)
par(mfrow=c(1,3))
for (sampleid in c("T1","T2","T3")){
  data=read.table(paste(sampleid,".combined.scaf.tsv",sep=""),header=T)
  #only look at sc1,2,3 variant
  data=data[c(1:2473),]
  if(sampleid=="T1"){
    colnames(data) = gsub(".1",".7",colnames(data))
    assign=read.table(paste(sampleid,".lymphocyte.assign",sep=""),sep='\t',header=T)
    assign = fillInQual(assign)
    assign$Barcode = gsub("-1",".7",assign$Barcode)
  }
  if(sampleid=="T2"){
    colnames(data) = gsub(".1",".8",colnames(data))
    assign=read.table(paste(sampleid,".lymphocyte.assign",sep=""),sep='\t',header=T)
    assign = fillInQual(assign)
    assign$Barcode = gsub("-1",".8",assign$Barcode)
  }
  if(sampleid=="T3"){
    colnames(data) = gsub(".1",".9",colnames(data))
    assign=read.table(paste(sampleid,".lymphocyte.assign",sep=""),sep='\t',header=T)
    assign = fillInQual(assign)
    assign$Barcode = gsub("-1",".9",assign$Barcode)
  }
  normal = assign[assign$ASIG=='normal',]$Barcode
  sc1 = assign[assign$ASIG=='SC1',]$Barcode
  sc2 = assign[assign$ASIG=='SC2',]$Barcode
  sc3 = assign[assign$ASIG=='SC3',]$Barcode
  unassign = assign[assign$ASIG=='UNASSIGN',]$Barcode
  is.all=T
  cells = idents_all[idents_all$cell.lineage=="lymphcytes",]$Barcode
  plotCellVariantBySubcloneSortByQualAndByCellType(cells=cells,celltype="lymphocytes",is.all=is.all,assign)
}
dev.off()

############ Plot genotype for CD34+ cells
ba_c11 = read.csv("../../../4pts_3tp_scRNAseq_analysis/Ba_C11_manual.csv")
ba_c11$Barcode = gsub("-",".",ba_c11$Barcode)
c11_cells = ba_c11[ba_c11$Ba_C11=="C11",]$Barcode
sampleid="T1"
data=read.table(paste(sampleid,".combined.scaf.tsv",sep=""),header=T)
data=data[c(1:2473),]
colnames(data) = gsub(".1",".7",colnames(data))
assign=read.table(paste(sampleid,".pre.assign",sep=""),sep='\t',header=T)
assign$Barcode = gsub("-1",".7",assign$Barcode)
assign=fillInQual(assign)
assign[assign$QUAL<7,]$ASIG = "UNASSIGN"

normal = assign[assign$ASIG=='normal',]$Barcode
sc1 = assign[assign$ASIG=='SC1',]$Barcode
sc2 = assign[assign$ASIG=='SC2',]$Barcode
sc3 = assign[assign$ASIG=='SC3',]$Barcode
unassign = assign[assign$ASIG=='UNASSIGN',]$Barcode
is.all=F
pdf("T1_CD34cells_subclone.pdf",width = 4, height=4)
plotCellVariantBySubcloneAndByCellType(cells=c11_cells,celltype="CD34+ cells",is.all=is.all)
dev.off()

umap=read.csv("../../../4pts_3tp_scRNAseq_analysis/4pts_3H_umapCoordinate.csv",header=T)
umap$Barcode = gsub("-",".",umap$Barcode)
rownames(umap) = umap$Barcode
umap = umap[rownames(umap) %in% c11_cells & 
              rownames(umap) %in% assign$Barcode,c(2,3)]
pdf("T1.CD34cells.subclones.umap.pdf",width=5,height=5)
plotUmap(xlim=c(-5,-2),ylim=c(2,3.5))
dev.off()