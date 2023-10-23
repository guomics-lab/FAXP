#In Tip vs In Gel####
peptideFiles<-list.files(pattern="peptide.tsv","1_InGelvsInTip/",full.names = T)
proteinFiles<-list.files(pattern="protein.tsv","1_InGelvsInTip/",full.names = T)

condFullNames<-unique(sapply(sapply(list.files("1_InGelvsInTip/"
                                               ,pattern = "tsv"),strsplit,"-"),"[[",1))

peptideMats<-list()
for(i in 1:length(peptideFiles)){
  peptideMats[[i]]<-read.delim(peptideFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}
proteinMats<-list()
for(i in 1:length(proteinFiles)){
  proteinMats[[i]]<-read.delim(proteinFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}
names(peptideMats)<-names(proteinMats)<-condFullNames

PepLen<-function(PepSeq){
  return(length(unlist(strsplit(PepSeq,""))))
}

pepList<-list()
for(i in 1:length(peptideMats)){
  pepList[[i]]<-(row.names(peptideMats[[i]]))
}
pepLen_InGel<-sapply(unique(unlist(pepList[1:6])),PepLen)
pepLen_InTip<-sapply(unique(unlist(pepList[7:12])),PepLen)

col2<-c("#FF9A00","#629AED")
names(col2)<-c("InTip","InGel")
pdf("1_InGelVsInTip_PepLen_line.pdf",height = 4,width=17)
PepLenFull_InTip<-rep(0,44)
PepLenFull_InGel<-rep(0,44)
names(PepLenFull_InTip)<-names(PepLenFull_InGel)<-seq(7,50)
PepLenFull_InTip[names(table(pepLen_InTip))]<-table(pepLen_InTip)
PepLenFull_InGel[names(table(pepLen_InGel))]<-table(pepLen_InGel)

# barplot(PepLenFull_InTip,col=col2["InTip"],space=0,border=NA,axes=F,ylim = c(0,4000))
# barplot(PepLenFull_InGel,add=T,col=col2["InGel"],space=0,border =NA,axes=F)

plot(PepLenFull_InTip,col=col2["InTip"],axes=F,ylim = c(0,4000),type='l'
     ,lwd=6)
lines(PepLenFull_InGel,add=T,col=col2["InGel"],lwd=6)

# hist(pepLen_InTip,col=col2["InTip"],border=NA,breaks = 100,space=0)
# hist(pepLen_InGel,add=T,col=col2["InGel"],border=NA)
# plot(as.numeric(names(table(pepLen_InTip))),table(pepLen_InTip),col=col2["InTip"],type='l',lwd=2.7,xlab="Peptide lengths",ylab="#",axes=F
     # ,ylim=c(0,4000))
# lines(as.numeric(names(table(pepLen_InGel))),table(pepLen_InGel),col=col2["InGel"],type='l',lwd=2.7)
axis(side=1,at=seq(7-6,50-6,3),labels = seq(7,50,3))
axis(side=2,at=seq(0,4000,1000),labels = seq(0,4000,1000),las=2)
legend("topright",legend=c("InTip","InGel"),fill=col2)
dev.off()

# This one is mainly used to visualize those sequences exclusively found under each conditions instead of all.
peptideFiles<-list.files(pattern="peptide.tsv","1_InGelvsInTip/",full.names = T)
proteinFiles<-list.files(pattern="protein.tsv","1_InGelvsInTip/",full.names = T)

condFullNames<-unique(sapply(sapply(list.files("1_InGelvsInTip/"
                                               ,pattern = "tsv"),strsplit,"-"),"[[",1))

peptideMats<-list()
for(i in 1:length(peptideFiles)){
  peptideMats[[i]]<-read.delim(peptideFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}
proteinMats<-list()
for(i in 1:length(proteinFiles)){
  proteinMats[[i]]<-read.delim(proteinFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}
names(peptideMats)<-names(proteinMats)<-condFullNames

PepLen<-function(PepSeq){
  return(length(unlist(strsplit(PepSeq,""))))
}

uniqueSeqsList <- list()
for (i in 1:length(peptideMats)) {
  uniqueSeqsList[[i]] <- unique(row.names(peptideMats[[i]]))
}
uniqueSeq_InGel <- unique(unlist(uniqueSeqsList[1:6]))
uniqueSeq_InTip <- unique(unlist(uniqueSeqsList[7:12]))
pepLen_Unique_InGel <- sapply(uniqueSeq_InGel, PepLen)
pepLen_Unique_InTip <- sapply(uniqueSeq_InTip, PepLen)
pepLen_OnlyInGel <- pepLen_Unique_InGel[!uniqueSeq_InGel %in% uniqueSeq_InTip]
pepLen_OnlyInTip <- pepLen_Unique_InTip[!uniqueSeq_InTip %in% uniqueSeq_InGel]

col2<-c("#FF9A00","#629AED")
names(col2)<-c("InTip","InGel")
pdf("1_InGelVsInTip_PepLen_line_dz.pdf",height = 4,width=4)
PepLenFull_InTip<-rep(0,44)
PepLenFull_InGel<-rep(0,44)
names(PepLenFull_InTip)<-names(PepLenFull_InGel)<-seq(7,50)
PepLenFull_InTip[names(table(pepLen_OnlyInTip))]<-table(pepLen_OnlyInTip)
PepLenFull_InGel[names(table(pepLen_OnlyInGel))]<-table(pepLen_OnlyInGel)

plot(PepLenFull_InGel,col=col2["InGel"],axes=F,ylim = c(0,4000),type='l'
     ,lwd=3)
lines(PepLenFull_InTip,add=T,col=col2["InTip"],lwd=3)
axis(side=1,at=seq(7-6,50-6,3),labels = seq(7,50,3))
axis(side=2,at=seq(0,2000,500),labels = seq(0,2000,500),las=2)
legend("topright",legend=c("InTip","InGel"),fill=col2)
dev.off()

df_InGelInTip<-as.data.frame(matrix(NA,nrow=length(condFullNames),ncol=3))
row.names(df_InGelInTip)<-condFullNames
colnames(df_InGelInTip)<-c("NumPeptide","NumProtein","MissClv")
df_InGelInTip$condName<-factor(sapply(sapply(condFullNames,strsplit,"_"),"[[",1))

df_InGelInTip$NumPeptide<-as.vector(sapply(peptideMats,nrow))
df_InGelInTip$NumProtein<-as.vector(sapply(proteinMats,nrow))

misClv<-function(PepSeq){
  if(length(unlist(gregexec("K|R",PepSeq)))>1){
    return(1)
  }else{
    return(0)
  }
}

pValueCal<-function(vectorIn,factorIn){
  levelNames<-levels(factorIn)
  comb_level = combn(levelNames,2)
  combNames<-paste(t(comb_level)[,1],t(comb_level)[,2],sep="-")
  pvalues<-c()
  for(j in 1:ncol(comb_level)){
    vector1<-vectorIn[factorIn==comb_level[1,j]]
    vector2<-vectorIn[factorIn==comb_level[2,j]]
    pvalues[j]<-t.test(vector1,vector2)$p.value
  }
  names(pvalues)<-combNames
  return(pvalues)
}

pValueCalSort<-function(vectorIn,factorIn){
  levelNames<-levels(factorIn)
  combNames<-paste(levelNames[1:(length(levelNames)-1)],levelNames[2:length(levelNames)])
  pvalues<-c()
  for(j in 1:(length(levelNames)-1)){
    vector1<-vectorIn[factorIn==levelNames[j]]
    vector2<-vectorIn[factorIn==levelNames[j+1]]
    pvalues[j]<-t.test(vector1,vector2)$p.value
  }
  names(pvalues)<-combNames
  return(pvalues)
}

pvalueStars<-function(pvalues){
  stars<-rep(NA,length(pvalues))
  stars[pvalues<=0.05&pvalues>0.01]<-"*"
  stars[pvalues<=0.01&pvalues>0.001]<-"**"
  stars[pvalues<=0.001&pvalues>0.0001]<-"***"
  stars[pvalues<=0.0001]<-"****"
  return(stars)
}


Mean_NumPeptide<-aggregate(df_InGelInTip$NumPeptide, list(df_InGelInTip$condName), FUN=mean)[,2]
Mean_NumProtein<-aggregate(df_InGelInTip$NumProtein, list(df_InGelInTip$condName), FUN=mean)[,2]

{pdf("1_InTipVsInGel_Id_colPvalue.pdf",width=5.2,height=3)
par(mfrow=c(1,2))
stripchart(df_InGelInTip$NumPeptide~df_InGelInTip$condName,col=col2,ylab="NumPeptide",vertical=T,method='jitter',pch=19,xlim=c(0.5,2.5)
           ,las=2,ylim=c(0,30000),cex=0)
segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_NumPeptide,y1=Mean_NumPeptide)
sem1<-aggregate(df_InGelInTip$NumPeptide, list(df_InGelInTip$condName), FUN=sd)[,2]
segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide+sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide-sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumPeptide+sem1,
         y1=Mean_NumPeptide+sem1)
stripchart(df_InGelInTip$NumPeptide~df_InGelInTip$condName,col=col2[levels(df_InGelInTip$condName)],ylab="NumPeptide",vertical=T,method='jitter',pch=19,las=2,add=T,cex=0.7)
pvalues<-pValueCal(df_InGelInTip$NumPeptide,df_InGelInTip$condName)
stars<-pvalueStars(pvalues)
arrows(x0=1,x1=2,y0=1.5*(Mean_NumPeptide[1]+Mean_NumPeptide[2])/2,y1=1.5*(Mean_NumPeptide[1]+Mean_NumPeptide[2])/2,code=3,angle=90,length = 0)
text(x=1.5,y=1.55*(Mean_NumPeptide[1]+Mean_NumPeptide[2])/2,labels=stars[!is.na(stars)])

stripchart(df_InGelInTip$NumProtein~df_InGelInTip$condName,col=col2,ylab="NumProtein",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,3500),cex=0,xlim=c(0.5,2.5))
segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_NumProtein,y1=Mean_NumProtein)
sem1<-aggregate(df_InGelInTip$NumProtein, list(df_InGelInTip$condName), FUN=sd)[,2]
segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_NumProtein-sem1,
         y1=Mean_NumProtein+sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumProtein-sem1,
         y1=Mean_NumProtein-sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumProtein+sem1,
         y1=Mean_NumProtein+sem1)
stripchart(df_InGelInTip$NumProtein~df_InGelInTip$condName,col=col2[levels(df_InGelInTip$condName)],ylab="NumProtein",vertical=T,method='jitter',pch=19,las=2,add=T,cex=0.7)
pvalues<-pValueCal(df_InGelInTip$NumProtein,df_InGelInTip$condName)
stars<-pvalueStars(pvalues)
arrows(x0=1,x1=2,y0=1.25*(Mean_NumProtein[1]+Mean_NumProtein[2])/2,y1=1.25*(Mean_NumProtein[1]+Mean_NumProtein[2])/2,code=3,angle=90,length = 0)
text(x=1.5,y=1.27*(Mean_NumProtein[1]+Mean_NumProtein[2])/2,labels=stars[!is.na(stars)])
dev.off()}

proList<-list()
for(i in 1:length(proteinMats)){
  proList[[i]]<-row.names(proteinMats[[i]])
}
Pep_InGel<-unique(unlist(pepList[1:6]))
Pep_InTip<-unique(unlist(pepList[7:12]))
Pro_InGel<-unique(unlist(proList[1:6]))
Pro_InTip<-unique(unlist(proList[7:12]))

InterPep<-intersect(Pep_InGel,Pep_InTip)
InterPro<-intersect(Pro_InGel,Pro_InTip)
Pep_bar<-c(length(Pep_InGel)-length(InterPep),length(InterPep),length(Pep_InTip)-length(InterPep))
Pro_bar<-c(length(Pro_InGel)-length(InterPro),length(InterPro),length(Pro_InTip)-length(InterPro))
pdf("1_InGelInTip_VennCircleCol.pdf",width = 4,height=4)
par(mfrow=c(3,1),mar=c(1,1,4,1))
grid::grid.newpage()
VennDiagram::draw.pairwise.venn(area1 = sum(Pep_bar[1:2]),area2 = sum(Pep_bar[2:3]),cross.area = Pep_bar[2],circle_order=c(T,T,T)
                                ,col =c(col2[2],col2[1]))
grid::grid.newpage()
VennDiagram::draw.pairwise.venn(area1 = sum(Pro_bar[1:2]),area2 = sum(Pro_bar[2:3]),cross.area = Pro_bar[2]
                                ,col =c(col2[2],col2[1]))
# bx<-barplot(as.matrix(Pep_bar),beside = F,horiz = T,col=c(col2[2],"#CBDB2A",col2[1]),main = "PepNum",axes=F)
# text(c(Pep_bar[1]/2,Pep_bar[1]+Pep_bar[2]/2,Pep_bar[1]+Pep_bar[2]+Pep_bar[3]/2),bx,(Pep_bar))
# 
# bx<-barplot(as.matrix(Pro_bar),beside = F,horiz = T,col=c(col2[2],"#CBDB2A",col2[1])
#         ,main = "ProNum",axes=F)
# text(c(Pro_bar[1]/2,Pro_bar[1]+Pro_bar[2]/2,Pro_bar[1]+Pro_bar[2]+Pro_bar[3]/2),bx,(Pro_bar))
# plot.new()
# legend("center",fill=c(col2[2],"#CBDB2A",col2[1]),legend = c("InGel Only","Intersective","InTip Only"))
dev.off()


#CV
Pep_InGel_Mat<-data.frame(matrix(NA,nrow=length(Pep_InGel),ncol=6))
Pep_InTip_Mat<-data.frame(matrix(NA,nrow=length(Pep_InTip),ncol=6))
Pro_InGel_Mat<-data.frame(matrix(NA,nrow=length(Pro_InGel),ncol=6))
Pro_InTip_Mat<-data.frame(matrix(NA,nrow=length(Pro_InTip),ncol=6))
colnames(Pep_InGel_Mat)<-colnames(Pro_InGel_Mat)<-condFullNames[1:6]
colnames(Pep_InTip_Mat)<-colnames(Pro_InTip_Mat)<-condFullNames[7:12]
row.names(Pep_InGel_Mat)<-Pep_InGel
row.names(Pep_InTip_Mat)<-Pep_InTip
row.names(Pro_InGel_Mat)<-Pro_InGel
row.names(Pro_InTip_Mat)<-Pro_InTip

for(i in 1:6){
  Pep_InGel_Mat[row.names(peptideMats[[i]]),i]<-peptideMats[[i]]$Intensity
  Pro_InGel_Mat[row.names(proteinMats[[i]]),i]<-proteinMats[[i]]$Total.Intensity
}
for(i in 7:12){
  Pep_InTip_Mat[row.names(peptideMats[[i]]),i-6]<-peptideMats[[i]]$Intensity
  Pro_InTip_Mat[row.names(proteinMats[[i]]),i-6]<-proteinMats[[i]]$Total.Intensity
}
Pep_InGel_CV<-apply(Pep_InGel_Mat,1,sd,na.rm=T)/apply(Pep_InGel_Mat,1,mean,na.rm=T)
Pep_InTip_CV<-apply(Pep_InTip_Mat,1,sd,na.rm=T)/apply(Pep_InTip_Mat,1,mean,na.rm=T)
Pro_InGel_CV<-apply(Pro_InGel_Mat,1,sd,na.rm=T)/apply(Pro_InGel_Mat,1,mean,na.rm=T)
Pro_InTip_CV<-apply(Pro_InTip_Mat,1,sd,na.rm=T)/apply(Pro_InTip_Mat,1,mean,na.rm=T)

CVList<-list(Pep_InGel_CV,Pep_InTip_CV,Pro_InGel_CV,Pro_InTip_CV)
names(CVList)<-c("Pep_InGel","Pep_InTip","Pro_InGel","Pro_InTip")
# pdf("1_InGelInTip_CV_col.pdf",width = 4,height = 4)
par(mar=c(6,4,1,1))
vioplot::vioplot(CVList,col=0,ylim=c(0,1),ylab="CoV",las=2)
boxplot(CVList,add=T,outline = F,col=col2[c("InGel","InTip")],las=2)
dev.off()

#gravy
gravyTable<-c(1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,3.8,-3.9,1.9,2.8,-1.6,-0.8,-0.7,-0.9,-1.3,4.2)
names(gravyTable)<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

gravyCal<-function(PepSeq){
  return(sum(gravyTable[unlist(strsplit(PepSeq,""))],na.rm=T))
}

Gravy_InGel<-sapply(row.names(Pep_InGel_Mat),gravyCal)
Gravy_InTip<-sapply(row.names(Pep_InTip_Mat),gravyCal)

pdf("1_InGelInTip_Gravy_col.pdf",width=6,height = 4)
plot(density(Gravy_InGel,bw=1.95),col=col2["InGel"],xlim=c(-40,30),xlab="GRAVY value",main="")
lines(density(Gravy_InTip,bw=1.95),col=col2["InTip"])
abline(v=median(Gravy_InGel),col=col2["InGel"])
abline(v=median(Gravy_InTip),col=col2["InTip"])
legend("topright",legend=c("InTip","InGel"),lwd=1.7,col=col2)
dev.off()

#cross correlation
pepAll<-unique(unlist(pepList))
proAll<-unique(unlist(proList))
Pep_Mat<-data.frame(matrix(NA,nrow=length(pepAll),ncol=12))
Pro_Mat<-data.frame(matrix(NA,nrow=length(proAll),ncol=12))
colnames(Pep_Mat)<-colnames(Pro_Mat)<-condFullNames
row.names(Pep_Mat)<-pepAll
row.names(Pro_Mat)<-proAll
for(i in 1:12){
  Pep_Mat[row.names(peptideMats[[i]]),i]<-peptideMats[[i]]$Intensity
  Pro_Mat[row.names(proteinMats[[i]]),i]<-proteinMats[[i]]$Total.Intensity
}
Pep_Mat_NoNA<-na.omit(Pep_Mat)
Pro_Mat_NoNA<-na.omit(Pro_Mat)
library(corrplot)
# CorList<-list(unique(as.vector(cor(preMatPoolNoNA)))[-1],unique(as.vector(cor(proMatPoolNoNA))[-1]))
PepCorMat<-(cor(Pep_Mat_NoNA))
ProCorMat<-(cor(Pro_Mat_NoNA))
PepCorMat_InGel<-PepCorMat[1:6,1:6]
ProCorMat_InGel<-ProCorMat[1:6,1:6]
PepCorMat_InTip<-PepCorMat[7:12,7:12]
ProCorMat_InTip<-ProCorMat[7:12,7:12]
pdf("Corr_Plot_col_0.pdf",width=6,height = 6)
toPlotList<-list(unique(as.vector(PepCorMat_InGel))[-1],
                 unique(as.vector(PepCorMat_InTip))[-1],
                 unique(as.vector(ProCorMat_InGel))[-1],
                 unique(as.vector(ProCorMat_InTip))[-1])
names(toPlotList)<-c("PepCor_InGel","PepCor_InTip",
                     "ProCor_InGel","ProCor_InTip")
boxplot(toPlotList,las=2,col=col2[c("InGel","InTip")],ylim=c(0,1))
dev.off()

# pdf("1_InGelInTip_Corr.pdf",width = 6,height = 6)
# par(mar=c(6,6,6,6))
# print(corrplot::corrplot(cor(Pep_Mat_NoNA),method = c("square"),tl.col=1,type="upper"
#          ,main="Peptides",addCoef.col="grey"))
# print(corrplot::corrplot(cor(Pro_Mat_NoNA),method = c("square"),tl.col=1,type="lower"
#          ,main="Proteins",add=T,addCoef.col="grey"))
# dev.off()

#mean abundance rank
ProteinInAll<-read.csv("191protein.csv",stringsAsFactors = F,header = T,row.names = 2)
ProteinIn<-ProteinInAll[,1]
AllProteinIn<-unique(sapply(sapply(unique(unlist(strsplit(ProteinIn,";"))),strsplit,"-"),"[[",1))
AllGeneIn<-c()
for(i in 1:length(AllProteinIn)){
  AllGeneIn[i]<-row.names(ProteinInAll)[which(grepl(AllProteinIn[i],ProteinInAll[,1]))]
  names(AllGeneIn)[i]<-AllProteinIn[i]
}
ProteinInUniq<-ProteinIn[!grepl(";",ProteinIn)]

pdf("1_InGelInTip_AbundRank_col_textGn.pdf",width = 8,height = 8)
par(mar=c(4,4,1,1))
plot(sort(log10(apply(Pro_InTip_Mat,1,mean,na.rm=T)),decreasing = T),type='l',col=col2["InTip"],lwd=2,xlim=c(0,4000)
     ,ylab="log10 mean protein intensities",xlab="Abundance rank",ylim=c(2,8))
ProSort<-sort(log10(apply(Pro_InTip_Mat,1,mean,na.rm=T)),decreasing = T)
ProSel<-which(sapply(sapply(names(ProSort),strsplit,"\\|"),"[[",2)%in%AllProteinIn)
points(ProSel,ProSort[ProSel],col=col2[1],pch=19)
text(ProSel+250,ProSort[ProSel],AllGeneIn[sapply(sapply(names(ProSort)[ProSel],strsplit,"\\|"),"[[",2)])
text(length(ProSort),min(ProSort[ProSort!="-Inf"]),round(min(ProSort[ProSort!="-Inf"]),2))
text(1,max(ProSort[ProSort!="-Inf"]),round(max(ProSort[ProSort!="-Inf"]),2))
lines(sort(log10(apply(Pro_InGel_Mat,1,mean,na.rm=T)),decreasing = T),type='l',col=col2["InGel"],lwd=2)
ProSort<-sort(log10(apply(Pro_InGel_Mat,1,mean,na.rm=T)),decreasing = T)
ProSel<-which(sapply(sapply(names(ProSort),strsplit,"\\|"),"[[",2)%in%AllProteinIn)
points(ProSel,ProSort[ProSel],col=col2[2],pch=19)
text(ProSel-250,ProSort[ProSel],AllGeneIn[sapply(sapply(names(ProSort)[ProSel],strsplit,"\\|"),"[[",2)])
text(length(ProSort),min(ProSort[ProSort!="-Inf"]),round(min(ProSort[ProSort!="-Inf"]),2))
text(1,max(ProSort[ProSort!="-Inf"]),round(max(ProSort[ProSort!="-Inf"]),2))
legend("topright",legend = c("InTip","InGel"),lwd=2,col = col2)
dev.off()

write.csv(row.names(Pro_InGel_Mat),"Pro_InGel.csv",row.names = F)
write.csv(row.names(Pro_InTip_Mat),"Pro_InTip.csv",row.names = F)

# Detection ranges####
#DIA
preMat_filename = '2_Detection_range_report.pr_matrix.tsv'
proMat_filename = '2_Detection_range_report.pg_matrix.tsv'

preMatIn<-read.table(preMat_filename,stringsAsFactors = F,sep='\t',header = T)
proMatIn<-read.table(proMat_filename,stringsAsFactors = F,sep='\t',header = T)

removeNames<-function(colnamein){
  return(paste(unlist(strsplit(unlist(strsplit(colnamein,"\\."))[12],"_"))[c(2:4,6)],collapse = "_"))
}
pepAll<-unique(preMatIn$Stripped.Sequence)
preCond<-sapply(colnames(preMatIn[,11:ncol(preMatIn)]),removeNames)
df_PepID<-as.data.frame(matrix(0,nrow=length(pepAll),ncol=length(preCond)))
colnames(df_PepID)<-preCond
row.names(df_PepID)<-pepAll

for(j in 11:ncol(preMatIn)){
  pepAllSel<-unique(preMatIn$Stripped.Sequence[!is.na(preMatIn[,j])])
  df_PepID[pepAllSel,j-10]<-1
}

proMat<-proMatIn[,6:ncol(proMatIn)]
row.names(proMat)<-proMatIn$Protein.Group
colnames(proMat)<-sapply(colnames(proMat),removeNames)
df_proID<-proMat
df_proID[!is.na(df_proID)]<-1
df_proID[is.na(df_proID)]<-0

DIA_pepNum<-apply(df_PepID,2,sum)
DIA_proNum<-apply(df_proID,2,sum)

DIA_Num<-as.data.frame(cbind(DIA_pepNum,DIA_proNum))

#DDA
#move files
FileLoc="Z:/members/Dongzhen/ProteomEx2_benchmark_result/N20230317ProteomEx_v2_detection_range_DDA"
dirAll<-list.dirs(FileLoc,full.names = F)[-1]
for(i in 1:length(dirAll)){
  file.copy(from=paste(FileLoc,dirAll[i],"protein.tsv",sep = "/"),
            to=paste0("2_Detection_range_DDA/",dirAll[i],"-protein.tsv"))
  file.copy(from=paste(FileLoc,dirAll[i],"peptide.tsv",sep = "/"),
            to=paste0("2_Detection_range_DDA/",dirAll[i],"-peptide.tsv"))
}

peptideFiles<-list.files(pattern="peptide.tsv","2_Detection_range_DDA/",full.names = T)
proteinFiles<-list.files(pattern="protein.tsv","2_Detection_range_DDA/",full.names = T)

peptideMats<-list()
for(i in 1:length(peptideFiles)){
  peptideMats[[i]]<-read.delim(peptideFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}
proteinMats<-list()
for(i in 1:length(proteinFiles)){
  proteinMats[[i]]<-read.delim(proteinFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}
names(peptideMats)<-names(proteinMats)<-condFullNames

pepList<-list()
for(i in 1:length(peptideMats)){
  pepList[[i]]<-(row.names(peptideMats[[i]]))
}
proList<-list()
for(i in 1:length(proteinMats)){
  proList[[i]]<-(row.names(proteinMats[[i]]))
}
condFullNames<-unique(sapply(sapply(list.files("2_Detection_range_DDA/"
                                               ,pattern = "tsv"),strsplit,"-"),"[[",1))
pepAll<-unique(unlist(pepList))
proAll<-unique(unlist(proList))
Pep_Mat<-data.frame(matrix(NA,nrow=length(pepAll),ncol=length(condFullNames)))
Pro_Mat<-data.frame(matrix(NA,nrow=length(proAll),ncol=length(condFullNames)))
colnames(Pep_Mat)<-colnames(Pro_Mat)<-condFullNames
row.names(Pep_Mat)<-pepAll
row.names(Pro_Mat)<-proAll
for(i in 1:length(condFullNames)){
  Pep_Mat[row.names(peptideMats[[i]]),i]<-peptideMats[[i]]$Intensity
  Pro_Mat[row.names(proteinMats[[i]]),i]<-proteinMats[[i]]$Total.Intensity
}
DDA_pepNum<-apply(!is.na(Pep_Mat),2,sum)
DDA_proNum<-apply(!is.na(Pro_Mat),2,sum)
DDA_Num<-as.data.frame(cbind(DDA_pepNum,DDA_proNum))
DDA_Num$CondName<-unlist(sapply(sapply(row.names(DDA_Num),strsplit,"_"),"[[",1))
DDA_Num$Diameter<-unlist(sapply(sapply(row.names(DDA_Num),strsplit,"_"),"[[",2))
DDA_Num$Diameter<-as.numeric(gsub("um","",gsub("mm","000um",DDA_Num$Diameter)))


DIA_Num$CondName<-unlist(sapply(sapply(row.names(DIA_Num),strsplit,"_"),"[[",1))
DIA_Num$Diameter<-unlist(sapply(sapply(row.names(DIA_Num),strsplit,"_"),"[[",2))
DIA_Num$Diameter<-as.numeric(gsub("um","",gsub("mm","000um",DIA_Num$Diameter)))

#
# stripchart(DIA_Num$DIA_proNum~DIA_Num$Diameter,ylab="proNum",vertical=T,method='jitter',pch=18
#            ,las=2,ylim=c(0,6000))
{
pdf("2_Detection_range_errCol_0.pdf",width=9,height = 5)
par(mfrow=c(1,2))
x000<-c(0.042,0.085,0.191,0.339,0.763,1.357,3.054)
conc<-c(350,500,750,1000,1500,2000,3000)

#peptide
DIA_Num_InTip<-DIA_Num[DIA_Num$CondName=="InTip",]
Mean<-aggregate(DIA_Num_InTip[,1], list(DIA_Num_InTip[,4]), FUN=mean)[,2]
Sd<-aggregate(DIA_Num_InTip[,1], list(DIA_Num_InTip[,4]), FUN=sd)[,2]
stripchart(DIA_Num_InTip[,1]~DIA_Num_InTip[,4],ylab="# peptides",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,50000),col=col2['InTip'],cex=0,at=x000,xlim=c(0,3000),axes=F)
segments(x0=x000-30,x1=x000+30,y0=Mean,y1=Mean,col=col2['InTip'])
segments(x0=x000,x1=x000,y0=Mean-Sd,y1=Mean+Sd,col=col2['InTip'])
segments(x0=x000-15,x1=x000+15,y0=Mean-Sd,y1=Mean-Sd,col=col2['InTip'])
segments(x0=x000-15,x1=x000+15,y0=Mean+Sd,y1=Mean+Sd,col=col2['InTip'])
# stripchart(DIA_Num_InTip[,1]~DIA_Num_InTip[,4],ylab="pepNum",vertical=T,method='jitter',pch=19
#            ,las=2,ylim=c(0,50000),col=col2['InTip'],add=T,cex=0.7,at=x000)

DIA_Num_InGel<-DIA_Num[DIA_Num$CondName=="InGel",]
Mean<-aggregate(DIA_Num_InGel[,1], list(DIA_Num_InGel[,4]), FUN=mean)[,2]
Sd<-aggregate(DIA_Num_InGel[,1], list(DIA_Num_InGel[,4]), FUN=sd)[,2]
x00<-x000[seq(4,7)]
segments(x0=x00-30,x1=x00+30,y0=Mean,y1=Mean,col=col2['InGel'])
segments(x0=x00,x1=x00,y0=Mean-Sd,y1=Mean+Sd,col=col2['InGel'])
segments(x0=x00-15,x1=x00+15,y0=Mean-Sd,y1=Mean-Sd,col=col2['InGel'])
segments(x0=x00-15,x1=x00+15,y0=Mean+Sd,y1=Mean+Sd,col=col2['InGel'])
# stripchart(DIA_Num_InGel[,1]~DIA_Num_InGel[,4],ylab="pepNum",vertical=T,method='jitter',pch=19
#            ,las=2,ylim=c(0,50000),col=col2['InGel'],add=T,cex=0.7,at=x00)

DDA_Num_InTip<-DDA_Num[DDA_Num$CondName=="InTip",]
Mean<-aggregate(DDA_Num_InTip[,1], list(DDA_Num_InTip[,4]), FUN=mean)[,2]
Sd<-aggregate(DDA_Num_InTip[,1], list(DDA_Num_InTip[,4]), FUN=sd)[,2]
x00<-x000[seq(1,7)]
segments(x0=x00-30,x1=x00+30,y0=Mean,y1=Mean,col=col2['InTip'])
segments(x0=x00,x1=x00,y0=Mean-Sd,y1=Mean+Sd,col=col2['InTip'])
segments(x0=x00-15,x1=x00+15,y0=Mean-Sd,y1=Mean-Sd,col=col2['InTip'])
segments(x0=x00-15,x1=x00+15,y0=Mean+Sd,y1=Mean+Sd,col=col2['InTip'])
# stripchart(DDA_Num_InTip[,1]~DDA_Num_InTip[,4],ylab="pepNum",vertical=T,method='jitter',pch=17
#            ,las=2,ylim=c(0,50000),col=col2['InTip'],add=T,cex=1,at=x00)


DDA_Num_InGel<-DDA_Num[DDA_Num$CondName=="InGel",]
Mean<-aggregate(DDA_Num_InGel[,1], list(DDA_Num_InGel[,4]), FUN=mean)[,2]
Sd<-aggregate(DDA_Num_InGel[,1], list(DDA_Num_InGel[,4]), FUN=sd)[,2]
x00<-x000[seq(4,7)]
segments(x0=x00-30,x1=x00+30,y0=Mean,y1=Mean,col=col2['InGel'])
segments(x0=x00,x1=x00,y0=Mean-Sd,y1=Mean+Sd,col=col2['InGel'])
segments(x0=x00-15,x1=x00+15,y0=Mean-Sd,y1=Mean-Sd,col=col2['InGel'])
segments(x0=x00-15,x1=x00+15,y0=Mean+Sd,y1=Mean+Sd,col=col2['InGel'])
# stripchart(DDA_Num_InGel[,1]~DDA_Num_InGel[,4],ylab="pepNum",vertical=T,method='jitter',pch=17
#            ,las=2,ylim=c(0,50000),col=col2['InGel'],add=T,cex=1,at=x00)
axis(side=3,at=c(0,x000),labels = c(0,x000),las=2)
axis(side=1,at=c(0,x000),labels = c(0,conc),las=2)
axis(side=2,at=seq(0,50000,10000),labels = seq(0,50000,10000),las=2)


lss=loess(c(0,DDA_Num_InGel[,1])~c(0,DDA_Num_InGel[,4]))
x00<-x000[seq(4,7)]
lines(c(0,seq(x00[1],x00[length(x00)],length.out=100)),c(0,predict(lss,seq(x00[1],x00[length(x00)],length.out=100)))
      ,col=col2['InGel'],lty=2)
lss=loess(c(0,DIA_Num_InGel[,1])~c(0,DIA_Num_InGel[,4]))
lines(c(0,seq(x00[1],x00[length(x00)],length.out=100)),c(0,predict(lss,seq(x00[1],x00[length(x00)],length.out=100)))
      ,col=col2['InGel'])
x00<-x000[seq(1,7)]
lss=loess(c(0,DDA_Num_InTip[,1])~c(0,DDA_Num_InTip[,4]))
lines(c(0,seq(x00[1],x00[length(x00)],length.out=100)),c(0,predict(lss,seq(x00[1],x00[length(x00)],length.out=100)))
      ,col=col2['InTip'],lty=2)
lss=loess(c(0,DIA_Num_InTip[,1])~c(0,DIA_Num_InTip[,4]))
lines(c(0,seq(x00[1],x00[length(x00)],length.out=100)),c(0,predict(lss,seq(x00[1],x00[length(x00)],length.out=100)))
      ,col=col2['InTip'])

#protein
DIA_Num_InTip<-DIA_Num[DIA_Num$CondName=="InTip",]
Mean<-aggregate(DIA_Num_InTip[,2], list(DIA_Num_InTip[,4]), FUN=mean)[,2]
Sd<-aggregate(DIA_Num_InTip[,2], list(DIA_Num_InTip[,4]), FUN=sd)[,2]
stripchart(DIA_Num_InTip[,2]~DIA_Num_InTip[,4],ylab="# proteins",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,6000),col=col2['InTip'],cex=0,at=x000,xlim=c(0,3000),axes=F)
segments(x0=x000-30,x1=x000+30,y0=Mean,y1=Mean,col=col2['InTip'])
segments(x0=x000,x1=x000,y0=Mean-Sd,y1=Mean+Sd,col=col2['InTip'])
segments(x0=x000-15,x1=x000+15,y0=Mean-Sd,y1=Mean-Sd,col=col2['InTip'])
segments(x0=x000-15,x1=x000+15,y0=Mean+Sd,y1=Mean+Sd,col=col2['InTip'])
# stripchart(DIA_Num_InTip[,2]~DIA_Num_InTip[,4],ylab="proNum",vertical=T,method='jitter',pch=19
#            ,las=2,ylim=c(3000,5500),col=col2['InTip'],add=T,cex=0.7,at=x000)

DIA_Num_InGel<-DIA_Num[DIA_Num$CondName=="InGel",]
Mean<-aggregate(DIA_Num_InGel[,2], list(DIA_Num_InGel[,4]), FUN=mean)[,2]
Sd<-aggregate(DIA_Num_InGel[,2], list(DIA_Num_InGel[,4]), FUN=sd)[,2]
x00<-x000[seq(4,7)]
segments(x0=x00-30,x1=x00+30,y0=Mean,y1=Mean,col=col2['InGel'])
segments(x0=x00,x1=x00,y0=Mean-Sd,y1=Mean+Sd,col=col2['InGel'])
segments(x0=x00-15,x1=x00+15,y0=Mean-Sd,y1=Mean-Sd,col=col2['InGel'])
segments(x0=x00-15,x1=x00+15,y0=Mean+Sd,y1=Mean+Sd,col=col2['InGel'])
# stripchart(DIA_Num_InGel[,2]~DIA_Num_InGel[,4],ylab="proNum",vertical=T,method='jitter',pch=19
#            ,las=2,ylim=c(0,50000),col=col2['InGel'],add=T,cex=0.7,at=x00)

DDA_Num_InTip<-DDA_Num[DDA_Num$CondName=="InTip",]
Mean<-aggregate(DDA_Num_InTip[,2], list(DDA_Num_InTip[,4]), FUN=mean)[,2]
Sd<-aggregate(DDA_Num_InTip[,2], list(DDA_Num_InTip[,4]), FUN=sd)[,2]
x00<-x000[seq(1,7)]
segments(x0=x00-30,x1=x00+30,y0=Mean,y1=Mean,col=col2['InTip'])
segments(x0=x00,x1=x00,y0=Mean-Sd,y1=Mean+Sd,col=col2['InTip'])
segments(x0=x00-15,x1=x00+15,y0=Mean-Sd,y1=Mean-Sd,col=col2['InTip'])
segments(x0=x00-15,x1=x00+15,y0=Mean+Sd,y1=Mean+Sd,col=col2['InTip'])
# stripchart(DDA_Num_InTip[,2]~DDA_Num_InTip[,4],ylab="proNum",vertical=T,method='jitter',pch=17
#            ,las=2,ylim=c(0,50000),col=col2['InTip'],add=T,cex=1,at=x00)


DDA_Num_InGel<-DDA_Num[DDA_Num$CondName=="InGel",]
Mean<-aggregate(DDA_Num_InGel[,2], list(DDA_Num_InGel[,4]), FUN=mean)[,2]
Sd<-aggregate(DDA_Num_InGel[,2], list(DDA_Num_InGel[,4]), FUN=sd)[,2]
x00<-x000[seq(4,7)]
segments(x0=x00-30,x1=x00+30,y0=Mean,y1=Mean,col=col2['InGel'])
segments(x0=x00,x1=x00,y0=Mean-Sd,y1=Mean+Sd,col=col2['InGel'])
segments(x0=x00-15,x1=x00+15,y0=Mean-Sd,y1=Mean-Sd,col=col2['InGel'])
segments(x0=x00-15,x1=x00+15,y0=Mean+Sd,y1=Mean+Sd,col=col2['InGel'])
# stripchart(DDA_Num_InGel[,2]~DDA_Num_InGel[,4],ylab="proNum",vertical=T,method='jitter',pch=17
#            ,las=2,ylim=c(0,50000),col=col2['InGel'],add=T,cex=1,at=x00)
axis(side=3,at=c(0,x000),labels = c(0,x000),las=2)
axis(side=1,at=c(0,x000),labels = c(0,conc),las=2)
axis(side=2,at=seq(0,6000,1000),labels = seq(0,6000,1000),las=2)

lss=loess(c(0,DDA_Num_InGel[,2])~c(0,DDA_Num_InGel[,4]))
x00<-x000[seq(4,7)]
lines(c(0,seq(x00[1],x00[length(x00)],length.out=100)),c(0,predict(lss,seq(x00[1],x00[length(x00)],length.out=100)))
      ,col=col2['InGel'],lty=2)
lss=loess(c(0,DIA_Num_InGel[,2])~c(0,DIA_Num_InGel[,4]))
lines(c(0,seq(x00[1],x00[length(x00)],length.out=100)),c(0,predict(lss,seq(x00[1],x00[length(x00)],length.out=100)))
      ,col=col2['InGel'])
x00<-x000[seq(1,7)]
lss=loess(c(0,DDA_Num_InTip[,2])~c(0,DDA_Num_InTip[,4]))
lines(c(0,seq(x00[1],x00[length(x00)],length.out=100)),c(0,predict(lss,seq(x00[1],x00[length(x00)],length.out=100)))
      ,col=col2['InTip'],lty=2)
lss=loess(c(0,DIA_Num_InTip[,2])~c(0,DIA_Num_InTip[,4]))
lines(c(0,seq(x00[1],x00[length(x00)],length.out=100)),c(0,predict(lss,seq(x00[1],x00[length(x00)],length.out=100)))
      ,col=col2['InTip'])
dev.off()
}

#minute sample####
preMat_filename = '3_minute_sample_report.pr_matrix.tsv'
proMat_filename = '3_minute_sample_report.pg_matrix.tsv'

preMatIn<-read.table(preMat_filename,stringsAsFactors = F,sep='/t',header = T)
proMatIn<-read.table(proMat_filename,stringsAsFactors = F,sep='/t',header = T)

removeNames<-function(colnamein){
  return(paste(unlist(strsplit(unlist(strsplit(colnamein,"//."))[12],"_"))[c(2:3,5)],collapse = "_"))
}
pepAll<-unique(preMatIn$Stripped.Sequence)
preCond<-sapply(colnames(preMatIn[,11:ncol(preMatIn)]),removeNames)
df_PepID<-as.data.frame(matrix(0,nrow=length(pepAll),ncol=length(preCond)))
colnames(df_PepID)<-preCond
row.names(df_PepID)<-pepAll

for(j in 11:ncol(preMatIn)){
  pepAllSel<-unique(preMatIn$Stripped.Sequence[!is.na(preMatIn[,j])])
  df_PepID[pepAllSel,j-10]<-1
}

proMat<-proMatIn[,6:ncol(proMatIn)]
row.names(proMat)<-proMatIn$Protein.Group
colnames(proMat)<-sapply(colnames(proMat),removeNames)
df_proID<-proMat
df_proID[!is.na(df_proID)]<-1
df_proID[is.na(df_proID)]<-0

Min_pepNum<-apply(df_PepID,2,sum)
Min_proNum<-apply(df_proID,2,sum)
Min_Num<-as.data.frame(cbind(Min_pepNum,Min_proNum))

pdf("3_minute_sample.pdf",width = 4,height = 4)
par(mfrow=c(1,2))
stripchart(Min_pepNum,ylab="pepNum",vertical=T,method='jitter',pch=17
           ,las=2,ylim=c(0,2000),cex=1)
Mean<-mean(Min_pepNum)
Sd<-sd(Min_pepNum)
x00<-seq(1,5)
segments(x0=x00-0.1,x1=x00+0.1,y0=Mean,y1=Mean)
segments(x0=x00,x1=x00,y0=Mean-Sd,y1=Mean+Sd)
segments(x0=x00-0.05,x1=x00+0.05,y0=Mean-Sd,y1=Mean-Sd)
segments(x0=x00-0.05,x1=x00+0.05,y0=Mean+Sd,y1=Mean+Sd)

stripchart(Min_proNum,ylab="proNum",vertical=T,method='jitter',pch=17
           ,las=2,ylim=c(0,600),cex=1)
Mean<-mean(Min_proNum)
Sd<-sd(Min_proNum)
x00<-seq(1,5)
segments(x0=x00-0.1,x1=x00+0.1,y0=Mean,y1=Mean)
segments(x0=x00,x1=x00,y0=Mean-Sd,y1=Mean+Sd)
segments(x0=x00-0.05,x1=x00+0.05,y0=Mean-Sd,y1=Mean-Sd)
segments(x0=x00-0.05,x1=x00+0.05,y0=Mean+Sd,y1=Mean+Sd)
dev.off()


#mis-cleavages
expNames<-list.files(pattern="exp","Z:/members/Dongzhen/ProteomEx2_demo_result/N20230322ProteomEx_v2_CRC_DDA/")
for(i in 1:5){
  file.copy(paste0("Z:/members/Dongzhen/ProteomEx2_demo_result/N20230322ProteomEx_v2_CRC_DDA/",expNames[i],"/peptide.tsv")
            ,paste0(expNames[i],"-peptide.tsv"))
  file.copy(paste0("Z:/members/Dongzhen/ProteomEx2_demo_result/N20230322ProteomEx_v2_CRC_DDA/",expNames[i],"/protein.tsv")
            ,paste0(expNames[i],"-protein.tsv"))
}

misClv<-function(PepSeq){
  if(length(unlist(gregexec("K|R",PepSeq)))>1){
    return(1)
  }else{
    return(0)
  }
}

misClvCRC<-c()
for(i in 1:5){
  readIn<-read.delim(paste0(expNames[i],"-peptide.tsv"),stringsAsFactors = F,sep='\t',row.names = 1,header = T)
  misClvCRC[i]<-mean(sapply(row.names(readIn),misClv))
}

for(i in 1:length(pepList)){
  df_InGelInTip$MissClv[i]<-mean(sapply(pepList[[i]], misClv))
}

pdf("MisClv_ML_CRC.pdf",width=7/3*2,height=4)
par(mfrow=c(1,2))
Mean_MissClv<-aggregate(df_InGelInTip$MissClv, list(df_InGelInTip$condName), FUN=mean)[,2]
stripchart(df_InGelInTip$MissClv~df_InGelInTip$condName,col=col2,ylab="MissClv",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0.2,0.5),cex=0,xlim=c(0.5,2.5),main="ML")
segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_MissClv,y1=Mean_MissClv)
sem1<-aggregate(df_InGelInTip$MissClv, list(df_InGelInTip$condName), FUN=sd)[,2]
segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_MissClv-sem1,
         y1=Mean_MissClv+sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MissClv-sem1,
         y1=Mean_MissClv-sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MissClv+sem1,
         y1=Mean_MissClv+sem1)
stripchart(df_InGelInTip$MissClv~df_InGelInTip$condName,col=col2[levels(df_InGelInTip$condName)],ylab="MissClv",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(5000,25000),add=T,cex=0.7)

stripchart(misClvCRC,col=1,ylab="MissClv",vertical=T,method='jitter',pch=19,las=2,cex=1,main="CRC")
Mean_MissClv<-mean(misClvCRC)
sem1<-sd(misClvCRC)
segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_MissClv,y1=Mean_MissClv)
segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_MissClv-sem1,
         y1=Mean_MissClv+sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MissClv-sem1,
         y1=Mean_MissClv-sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MissClv+sem1,
         y1=Mean_MissClv+sem1)
dev.off()

df_InGelInTip$MissClv<-c(20.40167364,21.58064263,21.27234855,21.45274453,34.4600401,23.64463099,17.54550094,17.10893227,17.53054833,16.88289947,16.9299836,17.46127108)
pdf("MisClv_ML_CRC_dz_colPvalue.pdf",width=7/3*2,height=4)
par(mfrow=c(1,2))
Mean_MissClv<-aggregate(df_InGelInTip$MissClv, list(df_InGelInTip$condName), FUN=mean)[,2]
stripchart(df_InGelInTip$MissClv~df_InGelInTip$condName,col=col2,ylab="Missed cleavage rate(%)",vertical=T,method='jitter',pch=19,ylim=c(10,40)
           ,las=2,cex=0,xlim=c(0.5,2.5),main="ML")
segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_MissClv,y1=Mean_MissClv)
sem1<-aggregate(df_InGelInTip$MissClv, list(df_InGelInTip$condName), FUN=sd)[,2]
segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_MissClv-sem1,
         y1=Mean_MissClv+sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MissClv-sem1,
         y1=Mean_MissClv-sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MissClv+sem1,
         y1=Mean_MissClv+sem1)
stripchart(df_InGelInTip$MissClv~df_InGelInTip$condName,col=col2[levels(df_InGelInTip$condName)],vertical=T,method='jitter',pch=19
           ,las=2,add=T,cex=0.7)
pvalues<-pValueCal(df_InGelInTip$MissClv,df_InGelInTip$condName)
stars<-pvalueStars(pvalues)
arrows(x0=1,x1=2,y0=1.8*(Mean_MissClv[1]+Mean_MissClv[2])/2,y1=1.8*(Mean_MissClv[1]+Mean_MissClv[2])/2,code=3,angle=90,length = 0)
text(x=1.5,y=1.81*(Mean_MissClv[1]+Mean_MissClv[2])/2,labels=stars[!is.na(stars)])

misClvCRC<-c(18.90,18.79,19.04,17.57,17.06)
stripchart(misClvCRC,col=1,ylab="Missed cleavage rate(%)",vertical=T,method='jitter',pch=19,las=2,cex=1,main="CRC",ylim=c(10,40))
Mean_MissClv<-mean(misClvCRC)
sem1<-sd(misClvCRC)
segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_MissClv,y1=Mean_MissClv)
segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_MissClv-sem1,
         y1=Mean_MissClv+sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MissClv-sem1,
         y1=Mean_MissClv-sem1)
segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MissClv+sem1,
         y1=Mean_MissClv+sem1)
dev.off()

#pathway enrichments
# InTip<-read.csv("InTip.csv",stringsAsFactors = F,row.names = 1,header = T)[2:20,]
# InGel<-read.csv("InGel.csv",stringsAsFactors = F,row.names = 1,header = T)[3:19,]
InTip_InGel_Loc<-t(read.csv("InTip_InGel_combinations_Location.csv",stringsAsFactors = F,row.names = 1,header = T))
InTip_InGel_Type<-t(read.csv("InTip_InGel_combinations_Type.csv",stringsAsFactors = F,row.names = 1,header = T))
row.names(InTip_InGel_Type)<-gsub("\\."," ",row.names(InTip_InGel_Type))
InTip_InGel_TypePerc<-InTip_InGel_Type
InTip_InGel_LocPerc<-InTip_InGel_Loc

InTip_InGel_TypePerc[,1]<-InTip_InGel_Type[,1]/sum(InTip_InGel_Type[,1])*100
InTip_InGel_TypePerc[,2]<-InTip_InGel_Type[,2]/sum(InTip_InGel_Type[,2])*100

InTip_InGel_LocPerc[,1]<-InTip_InGel_Loc[,1]/sum(InTip_InGel_Loc[,1])*100
InTip_InGel_LocPerc[,2]<-InTip_InGel_Loc[,2]/sum(InTip_InGel_Loc[,2])*100

pdf("PathEnrich.pdf",width=24,height=6)
par(mfrow=c(1,2),mar=c(5,6,1,2))
barplot(InTip_InGel_TypePerc,horiz = T,col=hcl.colors(nrow(InTip_InGel_Type), "Set 2"),las=2,main="InTip_InGel_Type",space=0,xlim=c(0,150))
legend("topright",legend = row.names(InTip_InGel_Type),fill=hcl.colors(nrow(InTip_InGel_Type), "Set 2"))
barplot(InTip_InGel_LocPerc,horiz = T,col=hcl.colors(nrow(InTip_InGel_Loc), "Set 3"),las=2,main="InTip_InGel_Loc",space=0
        ,xlim=c(0,150))
legend("topright",legend = row.names(InTip_InGel_Loc),fill=hcl.colors(nrow(InTip_InGel_Loc), "Set 3"))
dev.off()
