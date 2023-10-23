#batch design####
batchDesign<-read.csv("20230219_Batch_design_CRC_demo.csv",stringsAsFactors = F,row.names = 1)
SampleInfo<-read.csv("SampleInfo.csv",stringsAsFactors = F,header = F)

# batchDesign$NewID<-batchDesign$ID
# batchDesign$NewID[grepl("P1",batchDesign$NewID)]<-gsub("-","S",batchDesign$NewID[grepl("P1",batchDesign$NewID)])
# batchDesign$NewID[grepl("P2",batchDesign$NewID)]<-gsub("P2","P2S1",batchDesign$NewID[grepl("P2",batchDesign$NewID)])
# batchDesign$NewID[grepl("P3",batchDesign$NewID)]<-gsub("P3","P3S1",batchDesign$NewID[grepl("P3",batchDesign$NewID)])
# # table(batchDesign$Region[batchDesign$Patient=="P1"&batchDesign$Slide=="1"])
# 
# batchDesign1<-batchDesign[-1,]
# batchDesign$batchID<-sapply(sapply(row.names(batchDesign),strsplit,"_"),"[[",1)
# batchMatrix<-as.data.frame(matrix(NA,nrow=22,ncol=6),stringsAsFactors = F)
# colnames(batchMatrix)<-paste0("b",seq(1,6))
# row.names(batchMatrix)<-seq(1,22)
# batchMatrixPlot<-batchMatrix
# for(i in 1:nrow(batchDesign)){
#   batchID<-unlist(strsplit(row.names(batchDesign)[i],"_"))[1]
#   FileID<-unlist(strsplit(row.names(batchDesign)[i],"_"))[2]
#   batchMatrix[FileID,batchID]<-batchDesign$ID[i]
#   batchMatrixPlot[FileID,batchID]<-batchDesign$NewID[i]
# }
# # NewIDname=batchDesign$NewID[1]
# 
toPlotPar<-function(NewIDname){
  sss<-unlist(strsplit(NewIDname,""))
  patID<-paste0(sss[1:2],collapse = "")
  textID<-paste0(sss[3:length(sss)],collapse = "")
  subtypeID<-paste0(sss[5:(length(sss)-1)],collapse = "")
  return(c(patID,textID,subtypeID))
}
# 
# colSubtype<-c("grey","blue","purple","red2","red1","red3")
# names(colSubtype)<-c('N','L','H',"C","CC","PC")
# pchPat<-c(0,1,5)
# names(pchPat)<-c("P1","P2","P3")
# pdf("Vertical_batchMatrix.pdf",width=8,height = 12)
# plot(1,1,type="n",xlim=c(0,7),ylim=c(0,23),cex=2,lwd=2)
# for(i in 1:nrow(batchMatrixPlot)){
#   for(j in 1:ncol(batchMatrixPlot)){
#     plotPar<-toPlotPar(batchMatrixPlot[i,j])
#     lines(x=j,y=i,cex=5,type="p",lwd=2,
#     pch=pchPat[plotPar[1]],
#     col=colSubtype[plotPar[3]])
#     text(x=j,y=i,plotPar[2])
#   }
# }
# dev.off()
# pdf("Horizontal_batchMatrix.pdf",width=10,height = 5)
# plot(1,1,type="n",xlim=c(0,23),ylim=c(0,7),cex=2,lwd=2)
# for(i in 1:nrow(batchMatrixPlot)){
#   for(j in 1:ncol(batchMatrixPlot)){
#     plotPar<-toPlotPar(batchMatrixPlot[i,j])
#     lines(x=i,y=j,cex=3.2,type="p",lwd=2,
#           pch=pchPat[plotPar[1]],
#           col=colSubtype[plotPar[3]])
#     text(x=i,y=j,plotPar[2],cex=0.5)
#   }
# }
# dev.off()
# 
# pdf("Legend_batchMatrix.pdf")
# plot.new()
# legend("left",legend=names(colSubtype),col=colSubtype,lty=1)
# legend("center",legend=names(pchPat),pch=pchPat)
# dev.off()


#QC####
# preMat_filename = 'library_based_report.pr_matrix.tsv'
# proMat_filename = 'library_based_report.pg_matrix.tsv'
preMat_filename = '20230321Hybrid_lib_report.pr_matrix.tsv'
proMat_filename = '20230321Hybrid_lib_report.pg_matrix.tsv'
# preMat_filename = 'Whole_lib_report.pr_matrix.tsv'
# proMat_filename = 'Whole_lib_report.pg_matrix.tsv'

removeNames<-function(colnamein){
  return(paste(unlist(strsplit(unlist(strsplit(colnamein,"\\."))[7],"_"))[2:3],collapse = "_"))
  # return(paste(unlist(strsplit(unlist(strsplit(colnamein,"\\."))[13],"_"))[2:3],collapse = "_"))
}

preMatIn<-read.table(preMat_filename,stringsAsFactors = F,sep='\t',header = T)
proMatIn<-read.table(proMat_filename,stringsAsFactors = F,sep='\t',header = T)
preMat<-preMatIn[,11:ncol(preMatIn)]
row.names(preMat)<-preMatIn$Precursor.Id
colnames(preMat)<-sapply(colnames(preMat),removeNames)
proMat<-proMatIn[,6:ncol(proMatIn)]
row.names(proMat)<-proMatIn$Protein.Group
colnames(proMat)<-sapply(colnames(proMat),removeNames)


batchDesign<-read.csv("20230219_Batch_design_CRC_demo.csv",stringsAsFactors = F
                      ,row.names = 1)
batchDesign$NewID<-batchDesign$ID
batchDesign$NewID[grepl("P1",batchDesign$NewID)]<-gsub("-","S",batchDesign$NewID[grepl("P1",batchDesign$NewID)])
batchDesign$NewID[grepl("P2",batchDesign$NewID)]<-gsub("P2","P2S1",batchDesign$NewID[grepl("P2",batchDesign$NewID)])
batchDesign$NewID[grepl("P3",batchDesign$NewID)]<-gsub("P3","P3S1",batchDesign$NewID[grepl("P3",batchDesign$NewID)])
batchDesign1<-batchDesign[-1,]# first file missing
batchDesign1$batchID<-sapply(sapply(row.names(batchDesign1),strsplit,"_"),"[[",1)

preMatPool<-preMat[,grepl("pool",colnames(preMat))]
proMatPool<-proMat[,grepl("pool",colnames(proMat))]
preMatPoolNoNA<-na.omit(preMatPool)
proMatPoolNoNA<-na.omit(proMatPool)

preCV<-apply(preMatPoolNoNA,1,sd,na.rm=T)/apply(preMatPoolNoNA,1,mean,na.rm=T)
proCV<-apply(proMatPoolNoNA,1,sd,na.rm=T)/apply(proMatPoolNoNA,1,mean,na.rm=T)
pdf("QC_Pool_20230322.pdf")
par(mfrow=c(2,2),mar=c(4,4,6,4))

CVlist<-list(preCV,proCV)
names(CVlist)<-c("Precursors","Proteins")
vioplot::vioplot(CVlist,main="CV of pool samples",col=0)
boxplot(CVlist,main="CV of pool samples",col=0,add=T,outline = F)
library(corrplot)
CorList<-list(unique(as.vector(cor(preMatPoolNoNA)))[-1],unique(as.vector(cor(proMatPoolNoNA))[-1]))
names(CorList)<-c("Precursors","Proteins")
vioplot::vioplot(CorList,ylim=c(0.9,1),main="Pearson correlation of pool samples"
                 ,col=0)
boxplot(CorList,ylim=c(0.9,1),main="Pearson correlation of pool samples",add=T,col=0)

corrplot(cor(preMatPoolNoNA),method = "number",tl.col=1,type="upper"
         ,main="Precursor")
corrplot(cor(preMatPoolNoNA),method = "number",tl.col=1,type="lower"
         ,main="Proteins")
dev.off()

#Identifications####
preMatNoPool<-preMat[,row.names(batchDesign)[row.names(batchDesign)%in%colnames(preMat)]]
sum(is.na(preMatNoPool))/ncol(preMatNoPool)/nrow(preMatNoPool)

proMatNoPool<-proMat[,row.names(batchDesign)[row.names(batchDesign)%in%colnames(proMat)]]
sum(is.na(proMatNoPool))/ncol(proMatNoPool)/nrow(proMatNoPool)

batchDesign2<-batchDesign1
batchDesign2$Region<-factor(batchDesign2$Region
                            ,levels=c('N','L','H',"C","PC","CC"))
batchDesign2<-batchDesign2[order(batchDesign2$Patient, batchDesign2$Slide,batchDesign2$Region,batchDesign2$Rep),]

pdf("barplots_ID_20230322.pdf",width=8,height=12)
ProNum<-apply(!is.na(proMatNoPool),2,sum)
names(ProNum)<-paste(colnames(proMatNoPool),batchDesign1[colnames(proMatNoPool),"NewID"],sep = "-")
par(mar=c(12,6,3,1),mfrow=c(3,1))
barplot(sort(ProNum),las=2,cex.names = 0.37,space=0,col="#6667AB",ylim=c(0,6000))
abline(h=seq(1000,6000,1000),lty=2)
barplot(sort(ProNum),las=2,cex.names = 0.37,add=T,space=0,col="#6667AB",ylim=c(0,6000))

barplot((ProNum),las=2,cex.names = 0.37,space=0,col="#6667AB")
abline(h=seq(1000,6000,1000),lty=2)
barplot((ProNum),las=2,cex.names = 0.37,add=T,space=0,col="#6667AB")

proMatNoPool_sort<-proMatNoPool[,row.names(batchDesign2)[row.names(batchDesign2)%in%colnames(proMatNoPool)]]
ProNum_sort<-apply(!is.na(proMatNoPool_sort),2,sum)
names(ProNum_sort)<-paste(colnames(proMatNoPool_sort),batchDesign2[colnames(proMatNoPool_sort),"NewID"],sep = "-")
barplot((ProNum_sort),las=2,cex.names = 0.37,space=0,col="#6667AB")
abline(h=seq(1000,6000,1000),lty=2)
barplot((ProNum_sort),las=2,cex.names = 0.37,add=T,space=0,col="#6667AB")
dev.off()

#colorbars
# batchNewIDname=names(ProNum_sort)[1]
barPar<-function(batchNewIDname){
    batchNum<-unlist(strsplit(batchNewIDname,"_"))[1]
    plotPars<-toPlotPar(unlist(strsplit(batchNewIDname,"-"))[2])
    subType<-plotPars[3]
    PatID<-plotPars[1]
    SlideID<-paste(unlist(strsplit(plotPars[2],""))[1:2],collapse = "")
    return(c(batchNum,PatID,SlideID,subType))
}
# ProNum_sort
pdf("Colorbar_20230322.pdf",width=8,height=3)
GroupNum<-sapply(names(sort(ProNum)),barPar)
GroupNumFactor<-GroupNum
GroupNumFactor[1,]<-as.numeric(as.factor(GroupNum[1,]))
GroupNumFactor[2,]<-as.numeric(as.factor(GroupNum[2,]))*10+100
GroupNumFactor[3,]<-as.numeric(as.factor(GroupNum[3,]))*10+200
GroupNumFactor[4,]<-as.numeric(factor(GroupNum[4,],levels=c('N','L','H',"C","PC","CC")))*10+300
GroupNumMat<-matrix(0,nrow=nrow(GroupNumFactor),ncol=ncol(GroupNumFactor))
GroupNumMat[,]<-as.numeric(GroupNumFactor[,])
par(mfrow=c(4,1),mar=c(0,1,0,1))
image(t(t(GroupNumMat[4,])),col = hcl.colors(6, "viridis", rev = TRUE),axes=F)
image(t(t(GroupNumMat[3,])),col = hcl.colors(3, "Peach", rev = TRUE),axes=F)
image(t(t(GroupNumMat[2,])),col = hcl.colors(3, "Pastel 1", rev = TRUE),axes=F)
image(t(t(GroupNumMat[1,])),col = hcl.colors(6, "Spectral", rev = TRUE),axes=F)
dev.off()

# pdf("Legend_ID.pdf")
# plot.new()
# legend("topleft",legend=sort(unique(GroupNum[1,])),fill=hcl.colors(6, "Spectral", rev = TRUE))
# legend("bottomright",legend=sort(unique(GroupNum[2,])),fill=hcl.colors(3, "Pastel", rev = TRUE))
# legend("topright",legend=sort(unique(GroupNum[3,])),fill=hcl.colors(3, "Peach", rev = TRUE))
# legend("bottomleft",legend=sort(unique(batchDesign2$Region)),fill=hcl.colors(6, "viridis", rev = TRUE))
# dev.off()

BoxNumMat<-t(rbind(sort(ProNum),GroupNum))
colnames(BoxNumMat)<-c("ProNum","Batch","Patient","Slide","Subtype")
df_BoxNumMat<-as.data.frame(BoxNumMat)
df_BoxNumMat$ProNum<-as.numeric(df_BoxNumMat$ProNum)
df_BoxNumMat$Subtype<-factor(df_BoxNumMat$Subtype,levels=c('N','L','H',"C","PC","CC"))
pdf("Boxplots_ID_FZ20230322.pdf",width=10,height=3)
par(mfrow=c(1,4))
boxplot(df_BoxNumMat$ProNum~df_BoxNumMat$Batch,col=hcl.colors(6, "Spectral", rev = TRUE)
        ,xlab="",ylab="")
boxplot(df_BoxNumMat$ProNum~df_BoxNumMat$Patient,col=hcl.colors(3, "Pastel", rev = TRUE)
        ,xlab="",ylab="")
boxplot(df_BoxNumMat$ProNum~df_BoxNumMat$Slide,col=hcl.colors(3, "Peach", rev = TRUE)
        ,xlab="",ylab="")
boxplot(df_BoxNumMat$ProNum~df_BoxNumMat$Subtype,col=hcl.colors(6, "viridis", rev = TRUE)
        ,xlab="",ylab="")
dev.off()
# proMatPoolNA<-proMatPool[apply(is.na(proMatPool),1,sum)/ncol(proMatPool)<0.7,]
# uout0<-umap(t(proMatPoolNA))
# plot(uout0$layout[,1],uout0$layout[,2])

#overlaps####
proMatNoPool2<-proMatNoPool
colnames(proMatNoPool2)<-paste(colnames(proMatNoPool),batchDesign1[colnames(proMatNoPool),"NewID"],sep = "-")

proMatNoPool_P1<-proMatNoPool2[,grepl("P1",colnames(proMatNoPool2))]
proMatNoPool_P2<-proMatNoPool2[,grepl("P2",colnames(proMatNoPool2))]
proMatNoPool_P3<-proMatNoPool2[,grepl("P3",colnames(proMatNoPool2))]
proMatNoPool_P1<-proMatNoPool_P1[apply(is.na(proMatNoPool_P1),1,sum)!=ncol(proMatNoPool_P1),]
proMatNoPool_P2<-proMatNoPool_P2[apply(is.na(proMatNoPool_P2),1,sum)!=ncol(proMatNoPool_P2),]
proMatNoPool_P3<-proMatNoPool_P3[apply(is.na(proMatNoPool_P3),1,sum)!=ncol(proMatNoPool_P3),]

proMatNoPool_P1S1<-proMatNoPool2[,grepl("P1S1",colnames(proMatNoPool2))]
proMatNoPool_P1S2<-proMatNoPool2[,grepl("P1S2",colnames(proMatNoPool2))]
proMatNoPool_P1S3<-proMatNoPool2[,grepl("P1S3",colnames(proMatNoPool2))]
proMatNoPool_P1S1<-proMatNoPool_P1S1[apply(is.na(proMatNoPool_P1S1),1,sum)!=ncol(proMatNoPool_P1S1),]
proMatNoPool_P1S2<-proMatNoPool_P1S2[apply(is.na(proMatNoPool_P1S2),1,sum)!=ncol(proMatNoPool_P1S2),]
proMatNoPool_P1S3<-proMatNoPool_P1S3[apply(is.na(proMatNoPool_P1S3),1,sum)!=ncol(proMatNoPool_P1S3),]

# listIn<-list(row.names(proMatNoPool_P1),row.names(proMatNoPool_P2)
#              ,row.names(proMatNoPool_P3))
library(VennDiagram)
tripleVenn<-function(listIn,catNames){
  grid.newpage()
  a=listIn[[1]]
  b=listIn[[2]]
  c=listIn[[3]]
  draw.triple.venn(area1 = length(a),area2 = length(b),area3 = length(c),
                   n12 = length(intersect(a,b)),
                   n23 = length(intersect(b,c)),
                   n13 = length(intersect(a,c)),
                   n123 = length(intersect(a,intersect(b,c)))
                   ,fill = c("pink", "green", "orange"),
                   lty = "blank",
                   category = catNames)
}
pdf("Venn_P123-P1S123_20230322.pdf")
tripleVenn(list(row.names(proMatNoPool_P1),row.names(proMatNoPool_P2)
                             ,row.names(proMatNoPool_P3)),c("P1", "P2", "P3"))
tripleVenn(list(row.names(proMatNoPool_P1S1),row.names(proMatNoPool_P1S2)
                ,row.names(proMatNoPool_P1S3)),c("P1S1", "P1S2", "P1S3"))
dev.off()
library(UpSetR)
proMatNoPool_P1S1N<-proMatNoPool_P1S1[,colnames(proMatNoPool_P1S1)[grepl("N",colnames(proMatNoPool_P1S1))]]
proMatNoPool_P1S1L<-proMatNoPool_P1S1[,colnames(proMatNoPool_P1S1)[grepl("L",colnames(proMatNoPool_P1S1))]]
proMatNoPool_P1S1H<-proMatNoPool_P1S1[,colnames(proMatNoPool_P1S1)[grepl("H",colnames(proMatNoPool_P1S1))]]
proMatNoPool_P1S1PC<-proMatNoPool_P1S1[,colnames(proMatNoPool_P1S1)[grepl("PC",colnames(proMatNoPool_P1S1))]]
proMatNoPool_P1S1CC<-proMatNoPool_P1S1[,colnames(proMatNoPool_P1S1)[grepl("CC",colnames(proMatNoPool_P1S1))]]

proMatNoPool_P1S1N<-proMatNoPool_P1S1N[apply(is.na(proMatNoPool_P1S1N),1,sum)!=ncol(proMatNoPool_P1S1N),]
proMatNoPool_P1S1L<-proMatNoPool_P1S1L[apply(is.na(proMatNoPool_P1S1L),1,sum)!=ncol(proMatNoPool_P1S1L),]
proMatNoPool_P1S1H<-proMatNoPool_P1S1H[apply(is.na(proMatNoPool_P1S1H),1,sum)!=ncol(proMatNoPool_P1S1H),]
proMatNoPool_P1S1PC<-proMatNoPool_P1S1PC[apply(is.na(proMatNoPool_P1S1PC),1,sum)!=ncol(proMatNoPool_P1S1PC),]
proMatNoPool_P1S1CC<-proMatNoPool_P1S1CC[apply(is.na(proMatNoPool_P1S1CC),1,sum)!=ncol(proMatNoPool_P1S1CC),]

P1S1Pro<-row.names(proMatNoPool_P1S1)
P1S1ProMat<-proMatNoPool_P1S1[,1:5]
P1S1ProMat[,]<-0
colnames(P1S1ProMat)<-paste("P1S1",c("N","L","H","PC","CC"),sep = "_")
P1S1ProMat[row.names(proMatNoPool_P1S1N),1]<-1
P1S1ProMat[row.names(proMatNoPool_P1S1L),2]<-1
P1S1ProMat[row.names(proMatNoPool_P1S1H),3]<-1
P1S1ProMat[row.names(proMatNoPool_P1S1PC),4]<-1
P1S1ProMat[row.names(proMatNoPool_P1S1CC),5]<-1

library(UpSetR)
pdf("UpSet_P1S1_20230322.pdf",width= 8,height = 4)
upset(P1S1ProMat, decreasing = c(T,T),mainbar.y.label ="# Interescting proteins"
      ,sets.x.label = "# Proteins",show.numbers = "yes"
      ,order.by = c("freq"))
dev.off()
# upset(movies, nsets = 7, nintersects = NA, mb.ratio = c(0.5, 0.5),
     # order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
proMatNoPool_P1S1Mat<-proMatNoPool_P1S1
proMatNoPool_P1S1Mat[!is.na(proMatNoPool_P1S1)]<-1
proMatNoPool_P1S1Mat[is.na(proMatNoPool_P1S1)]<-0
# 
# pdf(paste0("UpSet_P1S1_20230321.pdf"),width= 8,height = 12)
# print(upset(proMatNoPool_P1S1Mat, decreasing = c(T,T),mainbar.y.label ="# Interescting proteins",sets.x.label = "# Proteins",nsets=50,show.numbers = "yes"
#             ,order.by = c("freq"),keep.order=T,mb.ratio=c(0.3,0.7)
#             ,nintersects = 30)
# )
# dev.off()

pdf(paste0("UpSetAll_20230322.pdf"),width= 8,height = 12)
proMatNoPool_Mat<-proMatNoPool
proMatNoPool_Mat[!is.na(proMatNoPool_Mat)]<-1
proMatNoPool_Mat[is.na(proMatNoPool_Mat)]<-0
colnames(proMatNoPool_Mat)<-paste(colnames(proMatNoPool),batchDesign1[colnames(proMatNoPool),"NewID"],sep = "-")
print(upset(proMatNoPool_Mat, decreasing = c(T,T),mainbar.y.label ="# Interescting proteins",sets.x.label = "# Proteins",nsets=131,show.numbers = "yes"
            ,order.by = c("freq"),keep.order=T,mb.ratio=c(0.3,0.7)
            ,text.scale = 0.6,nintersects = 50,point.size = 1.2,line.size = 0.3)
)
dev.off()

subTypes<-c('N','L','H',"PC","CC")
pdf(paste0("UpSet_P1S1_Subtypes_20230322.pdf"),width= 4,height = 4)
for(subType in subTypes){
  proMatNoPool_P1S1MatNoN<-proMatNoPool_P1S1Mat[,!grepl("N",colnames(proMatNoPool_P1S1Mat))]
  proMatNoPool_P1S1MatSub<-proMatNoPool_P1S1Mat[,grepl(subType,colnames(proMatNoPool_P1S1Mat))]
  # pdf(paste0("UpSet_P1S1_Subtypes",,"_20230313.pdf"),width= 4,height = 4)
  print(upset(proMatNoPool_P1S1MatSub, decreasing = c(T,T),mainbar.y.label ="# Interescting proteins",sets.x.label = "# Proteins",nsets=50,show.numbers = "yes"
        ,order.by = c("freq"),keep.order=T,mb.ratio=c(0.5,0.5)
        ,nintersects = 5)
)
}
dev.off()

write.csv(proMatNoPool,"proMatNoPool_FZ20230322.csv")

##Progression analysis of DEPs-added by dongzhen on 20231002####
batchDesign<-read.csv("20230219_Batch_design_CRC_demo.csv",stringsAsFactors = F,row.names = 1)
SampleInfo<-read.csv("SampleInfo.csv",stringsAsFactors = F,header = F)
batchDesign$NewID<-batchDesign$ID
batchDesign$NewID[grepl("P1",batchDesign$NewID)]<-gsub("-","S",batchDesign$NewID[grepl("P1",batchDesign$NewID)])
batchDesign$NewID[grepl("P2",batchDesign$NewID)]<-gsub("P2","P2S1",batchDesign$NewID[grepl("P2",batchDesign$NewID)])
batchDesign$NewID[grepl("P3",batchDesign$NewID)]<-gsub("P3","P3S1",batchDesign$NewID[grepl("P3",batchDesign$NewID)])
batchDesign1<-batchDesign[-1,]
batchDesign$batchID<-sapply(sapply(row.names(batchDesign),strsplit,"_"),"[[",1)
batchDesign2<-batchDesign1
batchDesign2$Region<-factor(batchDesign2$Region
                            ,levels=c('N','L','H',"C","PC","CC"))
batchDesign2<-batchDesign2[order(batchDesign2$Patient, batchDesign2$Slide,batchDesign2$Region,batchDesign2$Rep),]

subTypes<-c('N','L','H',"PC","CC")

proMatNoPool<-read.csv("proMatNoPool_DZ20231002.csv",stringsAsFactors = F,row.names = 1)
proMat<-proMatNoPool[!grepl(";",row.names(proMatNoPool)),]
colnames(proMat)<-batchDesign2[colnames(proMatNoPool),"NewID"]
quantile(apply(is.na(proMat),2,sum)/nrow(proMat))
quantile(apply(is.na(proMat),1,sum)/ncol(proMat))
quantile(apply(is.na(proMat),1,sum)/ncol(proMat),probs = seq(0, 1, 0.1))
#Choosing the threshold at a protein loss rate of around 0.5 can retain 70-80% of the protein.

##Progression analysis of DEPs of L vs. N-added by dongzhen on 20231002####
proMat_N<-proMat[,grepl("N",colnames(proMat))]
proMat_L<-proMat[,grepl("L",colnames(proMat))]
fc<-p<-c()
for(i in 1:nrow(proMat_L)){
  sel1<-as.numeric(proMat_L[i,])
  sel2<-as.numeric(proMat_N[i,])
  if(sum(!is.na(sel1))>2&sum(!is.na(sel2))>2){
    fc[i]<-mean(sel1,na.rm=T)/mean(sel2,na.rm=T)
    p[i]<-t.test(sel1,sel2)$p.value
  }else{
    fc[i]<- p[i]<-NA
  }
}
fc_LvsN<-fc
p_LvsN<-p
plot(log2(fc_LvsN),-log10(p_LvsN),xlim=c(-2,2))
abline(h=-log10(0.05),lty=2)
abline(v=log2(2^0.5),lty=2)
abline(v=log2(2^-0.5),lty=2)

p_LvsN_adjust<-p.adjust(p_LvsN,"BH")
up_LvsN<-row.names(proMat_L)[which(fc_LvsN>2^1&p_LvsN_adjust<0.05)]
dw_LvsN<-row.names(proMat_L)[which(fc_LvsN<2^-1&p_LvsN_adjust<0.05)]
names(fc_LvsN)<-names(p_LvsN)<-names(p_LvsN_adjust)<-row.names(proMat_L)

up_LvsN_filtered <- proMat_L[row.names(proMat_L) %in% up_LvsN, ]
dw_LvsN_filtered <- proMat_L[row.names(proMat_L) %in% dw_LvsN, ]

# Create a data frame containing filtered rows along with fc_LvsN and p_LvsN
up_LvsN_fc_p <- data.frame(RowName = row.names(up_LvsN_filtered), FoldChange = fc_LvsN[up_LvsN], PValue_adj = p_LvsN_adjust[up_LvsN])
dw_LvsN_fc_p <- data.frame(RowName = row.names(dw_LvsN_filtered), FoldChange = fc_LvsN[dw_LvsN], PValue_adj = p_LvsN_adjust[dw_LvsN])

# Save the combined data to LvsNup.csv with row names included
write.csv(up_LvsN_fc_p, "LvsNup.csv", row.names = FALSE)
write.csv(dw_LvsN_fc_p, "LvsNdw.csv", row.names = FALSE)

##Progression analysis of DEPs of H vs. N-added by dongzhen on 20231002####
proMat_N<-proMat[,grepl("N",colnames(proMat))]
proMat_H<-proMat[,grepl("H",colnames(proMat))]
fc<-p<-c()
for(i in 1:nrow(proMat_H)){
  sel1<-as.numeric(proMat_H[i,])
  sel2<-as.numeric(proMat_N[i,])
  if(sum(!is.na(sel1))>2&sum(!is.na(sel2))>2){
    fc[i]<-mean(sel1,na.rm=T)/mean(sel2,na.rm=T)
    p[i]<-t.test(sel1,sel2)$p.value
  }else{
    fc[i]<- p[i]<-NA
  }
}
fc_HvsN<-fc
p_HvsN<-p
plot(log2(fc_HvsN),-log10(p_HvsN),xlim=c(-2,2))
abline(h=-log10(0.05),lty=2)
abline(v=log2(2^0.5),lty=2)
abline(v=log2(2^-0.5),lty=2)

p_HvsN_adjust<-p.adjust(p_HvsN,"BH")
up_HvsN<-row.names(proMat_H)[which(fc_HvsN>2^1&p_HvsN_adjust<0.05)]
dw_HvsN<-row.names(proMat_H)[which(fc_HvsN<2^-1&p_HvsN_adjust<0.05)]
names(fc_HvsN)<-names(p_HvsN)<-names(p_HvsN_adjust)<-row.names(proMat_H)

up_HvsN_filtered <- proMat_H[row.names(proMat_H) %in% up_HvsN, ]
dw_HvsN_filtered <- proMat_H[row.names(proMat_H) %in% dw_HvsN, ]

# Create a data frame containing filtered rows along with fc_HvsN and p_HvsN
up_HvsN_fc_p <- data.frame(RowName = row.names(up_HvsN_filtered), FoldChange = fc_HvsN[up_HvsN], PValue_adj = p_HvsN_adjust[up_HvsN])
dw_HvsN_fc_p <- data.frame(RowName = row.names(dw_HvsN_filtered), FoldChange = fc_HvsN[dw_HvsN], PValue_adj = p_HvsN_adjust[dw_HvsN])

# Save the combined data to HvsNup.csv with row names included
write.csv(up_HvsN_fc_p, "HvsNup.csv", row.names = FALSE)
write.csv(dw_HvsN_fc_p, "HvsNdw.csv", row.names = FALSE)

##Progression analysis of DEPs of C vs. N-added by dongzhen on 20231002####
proMat_N<-proMat[,grepl("N",colnames(proMat))]
proMat_C<-proMat[,grepl("C",colnames(proMat))]
fc<-p<-c()
for(i in 1:nrow(proMat_C)){
  sel1<-as.numeric(proMat_C[i,])
  sel2<-as.numeric(proMat_N[i,])
  if(sum(!is.na(sel1))>2&sum(!is.na(sel2))>2){
    fc[i]<-mean(sel1,na.rm=T)/mean(sel2,na.rm=T)
    p[i]<-t.test(sel1,sel2)$p.value
  }else{
    fc[i]<- p[i]<-NA
  }
}
fc_CvsN<-fc
p_CvsN<-p
plot(log2(fc_CvsN),-log10(p_CvsN),xlim=c(-2,2))
abline(h=-log10(0.05),lty=2)
abline(v=log2(2^0.5),lty=2)
abline(v=log2(2^-0.5),lty=2)

p_CvsN_adjust<-p.adjust(p_HvsN,"BH")
up_CvsN<-row.names(proMat_C)[which(fc_CvsN>2^1&p_CvsN_adjust<0.05)]
dw_CvsN<-row.names(proMat_C)[which(fc_CvsN<2^-1&p_CvsN_adjust<0.05)]
names(fc_CvsN)<-names(p_CvsN)<-names(p_CvsN_adjust)<-row.names(proMat_C)

up_CvsN_filtered <- proMat_C[row.names(proMat_C) %in% up_CvsN, ]
dw_CvsN_filtered <- proMat_C[row.names(proMat_C) %in% dw_CvsN, ]

# Create a data frame containing filtered rows along with fc_CvsN and p_CvsN
up_CvsN_fc_p <- data.frame(RowName = row.names(up_CvsN_filtered), FoldChange = fc_CvsN[up_CvsN], PValue_adj = p_CvsN_adjust[up_CvsN])
dw_CvsN_fc_p <- data.frame(RowName = row.names(dw_CvsN_filtered), FoldChange = fc_CvsN[dw_CvsN], PValue_adj = p_CvsN_adjust[dw_CvsN])

# Save the combined data to CvsNup.csv with row names included
write.csv(up_CvsN_fc_p, "CvsNup.csv", row.names = FALSE)
write.csv(dw_CvsN_fc_p, "CvsNdw.csv", row.names = FALSE)

##Progression analysis of DEPs of H vs. L-added by dongzhen on 20231002####
proMat_L<-proMat[,grepl("L",colnames(proMat))]
proMat_H<-proMat[,grepl("H",colnames(proMat))]
fc<-p<-c()
for(i in 1:nrow(proMat_H)){
  sel1<-as.numeric(proMat_H[i,])
  sel2<-as.numeric(proMat_L[i,])
  if(sum(!is.na(sel1))>2&sum(!is.na(sel2))>2){
    fc[i]<-mean(sel1,na.rm=T)/mean(sel2,na.rm=T)
    p[i]<-t.test(sel1,sel2)$p.value
  }else{
    fc[i]<- p[i]<-NA
  }
}
fc_HvsL<-fc
p_HvsL<-p
plot(log2(fc_HvsL),-log10(p_HvsL),xlim=c(-2,2))
abline(h=-log10(0.05),lty=2)
abline(v=log2(2^0.5),lty=2)
abline(v=log2(2^-0.5),lty=2)

p_HvsL_adjust<-p.adjust(p_HvsL,"BH")
up_HvsL<-row.names(proMat_H)[which(fc_HvsL>2^1&p_HvsL_adjust<0.05)]
dw_HvsL<-row.names(proMat_H)[which(fc_HvsL<2^-1&p_HvsL_adjust<0.05)]
names(fc_HvsL)<-names(p_HvsL)<-names(p_HvsL_adjust)<-row.names(proMat_H)

up_HvsL_filtered <- proMat_H[row.names(proMat_H) %in% up_HvsL, ]
dw_HvsL_filtered <- proMat_H[row.names(proMat_H) %in% dw_HvsL, ]

# Create a data frame containing filtered rows along with fc_HvsL and p_HvsL
up_HvsL_fc_p <- data.frame(RowName = row.names(up_HvsL_filtered), FoldChange = fc_HvsL[up_HvsL], PValue_adj = p_HvsL_adjust[up_HvsL])
dw_HvsL_fc_p <- data.frame(RowName = row.names(dw_HvsL_filtered), FoldChange = fc_HvsL[dw_HvsL], PValue_adj = p_HvsL_adjust[dw_HvsL])

# Save the combined data to HvsLup.csv with row names included
write.csv(up_HvsL_fc_p, "HvsLup.csv", row.names = FALSE)
write.csv(dw_HvsL_fc_p, "HvsLdw.csv", row.names = FALSE)

##Progression analysis of DEPs of C vs. H-added by dongzhen on 20231002####
proMat_C<-proMat[,grepl("C",colnames(proMat))]
proMat_H<-proMat[,grepl("H",colnames(proMat))]
fc<-p<-c()
for(i in 1:nrow(proMat_C)){
  sel1<-as.numeric(proMat_C[i,])
  sel2<-as.numeric(proMat_H[i,])
  if(sum(!is.na(sel1))>2&sum(!is.na(sel2))>2){
    fc[i]<-mean(sel1,na.rm=T)/mean(sel2,na.rm=T)
    p[i]<-t.test(sel1,sel2)$p.value
  }else{
    fc[i]<- p[i]<-NA
  }
}
fc_CvsH<-fc
p_CvsH<-p
plot(log2(fc_CvsH),-log10(p_CvsH),xlim=c(-2,2))
abline(h=-log10(0.05),lty=2)
abline(v=log2(2^0.5),lty=2)
abline(v=log2(2^-0.5),lty=2)

p_CvsH_adjust<-p.adjust(p_CvsH,"BH")
up_CvsH<-row.names(proMat_C)[which(fc_CvsH>2^1&p_CvsH_adjust<0.05)]
dw_CvsH<-row.names(proMat_C)[which(fc_CvsH<2^-1&p_CvsH_adjust<0.05)]
names(fc_CvsH)<-names(p_CvsH)<-names(p_CvsH_adjust)<-row.names(proMat_C)

up_CvsH_filtered <- proMat_C[row.names(proMat_C) %in% up_CvsH, ]
dw_CvsH_filtered <- proMat_C[row.names(proMat_C) %in% dw_CvsH, ]

# Create a data frame containing filtered rows along with fc_CvsH and p_CvsH
up_CvsH_fc_p <- data.frame(RowName = row.names(up_CvsH_filtered), FoldChange = fc_CvsH[up_CvsH], PValue_adj = p_CvsH_adjust[up_CvsH])
dw_CvsH_fc_p <- data.frame(RowName = row.names(dw_CvsH_filtered), FoldChange = fc_CvsH[dw_CvsH], PValue_adj = p_CvsH_adjust[dw_CvsH])

# Save the combined data to CvsHup.csv with row names included
write.csv(up_CvsH_fc_p, "CvsHup.csv", row.names = FALSE)
write.csv(dw_CvsH_fc_p, "CvsHdw.csv", row.names = FALSE)

##Progression analysis of DEPs of CC vs. PC-added by dongzhen on 20231002####
proMat_P1_PC<-proMat[,grepl("PC",colnames(proMat))&grepl("P1",colnames(proMat))]
proMat_P1_CC<-proMat[,grepl("CC",colnames(proMat))&grepl("P1",colnames(proMat))]
fc<-p<-c()
for(i in 1:nrow(proMat_P1_CC)){
  sel1<-as.numeric(proMat_P1_CC[i,])
  sel2<-as.numeric(proMat_P1_PC[i,])
  if(sum(!is.na(sel1))>2&sum(!is.na(sel2))>2){
    fc[i]<-mean(sel1,na.rm=T)/mean(sel2,na.rm=T)
    p[i]<-t.test(sel1,sel2)$p.value
  }else{
    fc[i]<- p[i]<-NA
  }
}
fc_CCvsPC<-fc
p_CCvsPC<-p
plot(log2(fc_P1),-log10(p_P1),xlim=c(-2,2))
abline(h=-log10(0.05),lty=2)
abline(v=log2(2^0.5),lty=2)
abline(v=log2(2^-0.5),lty=2)

p_CCvsPC_adjust<-p.adjust(p_CCvsPC,"BH")
up_CCvsPC<-row.names(proMat_P1_CC)[which(fc_CCvsPC>2^1&p_CCvsPC_adjust<0.05)]
dw_CCvsPC<-row.names(proMat_P1_CC)[which(fc_CCvsPC<2^-1&p_CCvsPC_adjust<0.05)]
names(fc_CCvsPC)<-names(p_CCvsPC)<-names(p_CCvsPC_adjust)<-row.names(proMat_P1_CC)

up_CCvsPC_filtered <- proMat_P1_CC[row.names(proMat_P1_CC) %in% up_CCvsPC, ]
dw_CCvsPC_filtered <- proMat_P1_CC[row.names(proMat_P1_CC) %in% dw_CCvsPC, ]

# Create a data frame containing filtered rows along with fc_CCvsPC and p_CCvsPC
up_CCvsPC_fc_p <- data.frame(RowName = row.names(up_CCvsPC_filtered), FoldChange = fc_CCvsPC[up_CCvsPC], PValue_adj = p_CCvsPC_adjust[up_CCvsPC])
dw_CCvsPC_fc_p <- data.frame(RowName = row.names(dw_CCvsPC_filtered), FoldChange = fc_CCvsPC[dw_CCvsPC], PValue_adj = p_CCvsPC_adjust[dw_CCvsPC])

# Save the combined data to CCvsPCup.csv with row names included
write.csv(up_CCvsPC_fc_p, "CCvsPCup.csv", row.names = FALSE)
write.csv(dw_CCvsPC_fc_p, "CCvsPCdw.csv", row.names = FALSE)

## UpSet viaualization
library(UpSetR)
proMat_all<-read.csv("DEPs_up_dw20231003.csv",stringsAsFactors = F)

# Merge all rows in six columns and remove duplicate rows
combined_rows <- unique(unlist(proMat_all))

# Remove empty character lines
combined_rows <- combined_rows[combined_rows != ""]

# Create a new data frame with the merged rows as row names
result_df <- data.frame(matrix(0, nrow = length(combined_rows), ncol = ncol(proMat_all)))

# Set new row name
rownames(result_df) <- combined_rows

# Iterate over each column
for (col in colnames(proMat_all)) {
  # Get the unique element in the current column
  unique_elements <- unique(proMat_all[[col]])
  
  # Iterate through each row of the data frame
  for (i in 1:nrow(result_df)) {
    # Flag 1 if the string of the current row is in the list of unique elements of the current column
    if (rownames(result_df)[i] %in% unique_elements) {
      result_df[i, col] <- 1
    }
  }
}

result_df[is.na(result_df)] <- 0
result_df2 <- result_df[, (ncol(result_df) - 5):ncol(result_df)]
#write.csv(result_df2, "result_df2.csv", row.names = TRUE)

pdf("UpSet_DEPs_20231003.pdf",width= 8,height = 5)
upset(result_df2, decreasing = c(T,T),mainbar.y.label ="# Interescting proteins"
      ,sets.x.label = "# Proteins",show.numbers = "yes"
      ,nsets = 6
      ,order.by = c("freq"))
dev.off()