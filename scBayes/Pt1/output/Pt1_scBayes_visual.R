
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
plotCellVariantBySubcloneSortByQualAndByCellType = function(cells, celltype,is.all,assignment){
  order = order(as.double(assignment$QUAL),decreasing=T)
  ordered_data = data[,c(1,order+1)]
  print(celltype)
  sc1_cell_index = colnames(ordered_data) %in% cells & (colnames(ordered_data) %in% sc1)
  normal_index = colnames(ordered_data) %in% cells & (colnames(ordered_data) %in% normal)
  unassign_index = colnames(ordered_data) %in% cells & (colnames(ordered_data) %in% unassign)
  sc1_data = ordered_data[,sc1_cell_index]
  print("sc1 cell number")
  print(sum(sc1_cell_index))
  print("normal cells")
  print(sum(normal_index))
  print("unassigned cells")
  print(sum(unassign_index))
  normal_data = ordered_data[,normal_index]
  unassign_data = ordered_data[,unassign_index]
  
  if (is.all){
    new_data = cbind(ordered_data$variant,sc1_data,normal_data,unassign_data)
  } else {
    new_data = cbind(ordered_data$variant,sc1_data,normal_data)
  }
  colnames(new_data)[1]="variants"
  plotCellVariantMatrix(new_data,mainlabel=paste("genotype for ",celltype,sep=""),colnames(new_data)[2:ncol(new_data)])
  abline(v = sum(sc1_cell_index)+0.5)
  if (is.all) {
    abline(v = sum(sc1_cell_index)+0.5)
    abline(v = sum(sc1_cell_index)+ncol(normal_data)+0.5)  }
}
plotCellVariantBySubcloneAndByCellType = function(cells, celltype,is.all){
  print(celltype)
  tumor_cell_index = colnames(data) %in% cells & (colnames(data) %in% tumor)
  normal_index = colnames(data) %in% cells & (colnames(data) %in% normal)
  unassign_index = colnames(data) %in% cells & (colnames(data) %in% unassign)
  tumor_data = data[,tumor_cell_index]
  print("tumor cell number")
  print(sum(tumor_cell_index))
  print("normal cells")
  print(sum(normal_index))
  print("unassigned cells")
  print(sum(unassign_index))
  normal_data = data[,normal_index]
  unassign_data = data[,unassign_index]
  if (is.all){
    new_data = cbind(data$variant,tumor_data,normal_data,unassign_data)
  } else {
    new_data = cbind(data$variant,tumor_data)
  }
  colnames(new_data)[1]="variants"
  plotCellVariantMatrix(new_data,mainlabel=paste("genotype for ",celltype,sep=""),colnames(new_data)[2:ncol(new_data)])
  if (is.all) {
    abline(v = sum(tumor_cell_index)+0.5)
    abline(v = sum(tumor_cell_index)+ncol(normal_data)+0.5)
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

########## plot myeloid cells from three timepoint together
idents_all=read.csv("../../../4pts_3tp_scRNAseq_analysis/cell_lineage.csv",header=T)
idents_all$Barcode = gsub("-",".",idents_all$Barcode)

pdf("myeloid cells.pdf",width=13,height=4)
par(mfrow=c(1,3))

for (sampleid in c("T1","T2","T3")){
  data=read.table(paste(sampleid,".strictSomatic.scaf.tsv",sep=""),header=T)
  if(sampleid=="T1"){
    colnames(data) = gsub(".1",".10",colnames(data))
    assign=read.table(paste(sampleid,".pre.assign",sep=""),sep='\t',header=T)
    assign$Barcode = gsub("-1",".10",assign$Barcode)
    assign = fillInQual(assign)
    assign[assign$QUAL<14,]$ASIG = "UNASSIGN"
  }
  if(sampleid=="T2"){
    colnames(data) = gsub(".1",".11",colnames(data))
    assign=read.table(paste(sampleid,".post.assign",sep=""),sep='\t',header=T)
    assign$Barcode = gsub("-1",".11",assign$Barcode)
    assign = fillInQual(assign)
    assign[assign$QUAL<14,]$ASIG = "UNASSIGN"
    
  }
  if(sampleid=="T3"){
    colnames(data) = gsub(".1",".12",colnames(data))
    assign=read.table(paste(sampleid,".post.assign",sep=""),sep='\t',header=T)
    assign$Barcode = gsub("-1",".12",assign$Barcode)
    assign = fillInQual(assign)
    assign[assign$QUAL<14,]$ASIG = "UNASSIGN"
    
  }
  normal = assign[assign$ASIG=='normal',]$Barcode
  sc1 = assign[assign$ASIG=='tumor',]$Barcode
  unassign = assign[assign$ASIG=='UNASSIGN',]$Barcode
  is.all=T
  cells = idents_all[idents_all$cell.lineage=="myeloid_cells",]$Barcode
  plotCellVariantBySubcloneSortByQualAndByCellType(cells=cells,celltype="myeloid cells",is.all=is.all,assign)
}
dev.off()

########## plot lymphocytes cells from three timepoint together
pdf("lymphocytes.pdf",width=13,height=4)
par(mfrow=c(1,3))

for (sampleid in c("T1","T2","T3")){
  data=read.table(paste(sampleid,".strictSomatic.scaf.tsv",sep=""),header=T)
  assign=read.table(paste(sampleid,".lymphocyte.assign",sep=""),sep='\t',header=T)
  assign = fillInQual(assign)
  if(sampleid=="T1"){
    colnames(data) = gsub(".1",".10",colnames(data))
    assign=read.table(paste(sampleid,".lymphocyte.assign",sep=""),sep='\t',header=T)
    assign$Barcode = gsub("-1",".10",assign$Barcode)
  }
  if(sampleid=="T2"){
    colnames(data) = gsub(".1",".11",colnames(data))
    assign=read.table(paste(sampleid,".lymphocyte.assign",sep=""),sep='\t',header=T)
    assign$Barcode = gsub("-1",".11",assign$Barcode)
  }
  if(sampleid=="T3"){
    colnames(data) = gsub(".1",".12",colnames(data))
    assign=read.table(paste(sampleid,".lymphocyte.assign",sep=""),sep='\t',header=T)
    assign$Barcode = gsub("-1",".12",assign$Barcode)
  }
  normal = assign[assign$ASIG=='normal',]$Barcode
  sc1 = assign[assign$ASIG=='tumor',]$Barcode
  unassign = assign[assign$ASIG=='UNASSIGN',]$Barcode
  is.all=T
  cells = idents_all[idents_all$cell.lineage=="lymphcytes",]$Barcode
  plotCellVariantBySubcloneSortByQualAndByCellType(cells=cells,celltype="lymphocytes",is.all=is.all,assign)
}
dev.off()
