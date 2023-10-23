# The data is located in the Figure 1 folder under the "All_matrix" folder.
# unique(sapply(sapply(list.files("."),strsplit,"_"),"[[",2))[1:3]
# ReductionAlkylations####
peptideFiles<-list.files(pattern="peptide.tsv","1_ReductionAlkylations/",full.names = T)
proteinFiles<-list.files(pattern="protein.tsv","1_ReductionAlkylations/",full.names = T)

condFullNames<-unique(sapply(sapply(list.files("1_ReductionAlkylations/"
                                               ,pattern = "tsv"),strsplit,"-"),"[[",1))
condNames<-unique(sapply(sapply(condFullNames,strsplit,"_"),"[[",1))
names(condFullNames)<-sapply(sapply(condFullNames,strsplit,"_"),"[[",1)
condCol<-hcl.colors(length(condNames), palette = "Dynamic", alpha = 0.9)
names(condCol)<-condNames
# "#DB9D85E6" "#ABB065E6" "#5CBD92E6" "#4CB9CCE6" "#ACA4E2E6" "#E093C3E6"
peptideMats<-list()
for(i in 1:length(peptideFiles)){
  peptideMats[[i]]<-read.delim(peptideFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}

proteinMats<-list()
for(i in 1:length(proteinFiles)){
  proteinMats[[i]]<-read.delim(proteinFiles[i],stringsAsFactors = F,header = T,sep='\t',row.names = 1)
}
names(peptideMats)<-names(proteinMats)<-condFullNames

df_ReducAkyl<-as.data.frame(matrix(NA,nrow=length(condFullNames),ncol=4))
row.names(df_ReducAkyl)<-condFullNames
colnames(df_ReducAkyl)<-c("NumPeptide","NumProtein","CysModPerc","CysPeptidePerc")
df_ReducAkyl$condName<-names(condFullNames)
df_ReducAkyl$Rep<-c(rep(c(1,1,2,2,3,3),4),rep(1,7))

df_ReducAkyl$NumPeptide<-as.vector(sapply(peptideMats,nrow))
df_ReducAkyl$NumProtein<-as.vector(sapply(proteinMats,nrow))

cysMod<-function(pepMat){
  cysModNum<-length(which(grepl("C\\(57.0214\\)",pepMat$Assigned.Modifications)))
  cysPepNum<-sum(grepl("C",row.names(pepMat)))
  return(list(round(cysModNum/cysPepNum*100,2),round(cysPepNum/nrow(pepMat)*100,2)))
}
df_ReducAkyl$CysModPerc<-unlist(sapply(peptideMats,cysMod)[1,])
df_ReducAkyl$CysPeptidePerc<-unlist(sapply(peptideMats,cysMod)[2,])
write.csv(df_ReducAkyl,"1_ReductionAlkylations.csv")

df_ReducAkyl$condName<-as.factor(df_ReducAkyl$condName)
par(mfrow=c(1,4),mar=c(8,4,2,2))
stripchart(NumPeptide~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,34000))
stripchart(NumProtein~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,3400))
stripchart(CysModPerc~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(50,100))
stripchart(CysPeptidePerc~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,15))

#noMLQC
df_ReducAkylNoML<-df_ReducAkyl[!grepl("QC",row.names(df_ReducAkyl)),]
df_ReducAkylNoML$condName<-as.vector(df_ReducAkylNoML$condName)
df_ReducAkylNoML$condName<-as.factor(df_ReducAkylNoML$condName)
df_ReducAkylNoML$condName<-factor(df_ReducAkylNoML$condName,
                                  levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1"))

Mean_NumPeptide<-aggregate(df_ReducAkylNoML$NumPeptide, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_NumProtein<-aggregate(df_ReducAkylNoML$NumProtein, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_CysModPerc<-aggregate(df_ReducAkylNoML$CysModPerc, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_CysPeptidePerc<-aggregate(df_ReducAkylNoML$CysPeptidePerc, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
names(Mean_NumPeptide)<-names(Mean_NumProtein)<-names(Mean_CysModPerc)<-names(Mean_CysPeptidePerc)<-c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")



df_ReducAkylNoML4<-df_ReducAkylNoML3<-df_ReducAkylNoML2<-df_ReducAkylNoML1<-df_ReducAkylNoML
df_ReducAkylNoML1$condName<-factor(df_ReducAkylNoML1$condName,levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")[order(Mean_NumPeptide,decreasing = T)])
df_ReducAkylNoML2$condName<-factor(df_ReducAkylNoML2$condName,levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")[order(Mean_NumProtein,decreasing = T)])
df_ReducAkylNoML3$condName<-factor(df_ReducAkylNoML3$condName,levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")[order(Mean_CysModPerc,decreasing = T)])
df_ReducAkylNoML4$condName<-factor(df_ReducAkylNoML4$condName,levels=c("20mMTCEP","45mMTCEP","10mMDTT","20mMDTT","InGelv1")[order(Mean_CysPeptidePerc,decreasing = T)])

pValueCal<-function(vectorIn,factorIn){
  # vectorIn=df_ReducAkylNoML3$NumPeptide
  # factorIn=df_ReducAkylNoML3$condName
  
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
  # vectorIn=df_ReducAkylNoML3$CysModPerc
  # factorIn=df_ReducAkylNoML3$condName
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

{
pdf("1_ReductionAlkylations_pvalue.pdf",width=8.27,height=2.71)
par(mfrow=c(1,4),mar=c(8,4,1,1))
stripchart(df_ReducAkylNoML3$CysModPerc[df_ReducAkylNoML3$Rep==1]~df_ReducAkylNoML3$condName[df_ReducAkylNoML3$Rep==1],vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(94,100),col=condCol[levels(df_ReducAkylNoML3$condName)],ylab="CysModPerc")
segments(x0=seq(1,5)-0.1,x1=seq(1,5)+0.1,y0=sort(Mean_CysModPerc,decreasing = T),y1=sort(Mean_CysModPerc,decreasing = T))
sem3<-aggregate(df_ReducAkylNoML[,3], list(df_ReducAkylNoML3$condName), FUN=sd)[,2]
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_CysModPerc,decreasing = T)-sem3,
         y1=sort(Mean_CysModPerc,decreasing = T)+sem3)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_CysModPerc,decreasing = T)-sem3,
         y1=sort(Mean_CysModPerc,decreasing = T)-sem3)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_CysModPerc,decreasing = T)+sem3,
         y1=sort(Mean_CysModPerc,decreasing = T)+sem3)
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_CysModPerc,decreasing = T)-sem3,
         y1=sort(Mean_CysModPerc,decreasing = T)+sem3)
stripchart(df_ReducAkylNoML3$CysModPerc[df_ReducAkylNoML3$Rep==1]~df_ReducAkylNoML3$condName[df_ReducAkylNoML3$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(94,100),col=condCol[levels(df_ReducAkylNoML3$condName)],ylab="CysModPerc")
stripchart(df_ReducAkylNoML3$CysModPerc[df_ReducAkylNoML3$Rep==2]~df_ReducAkylNoML3$condName[df_ReducAkylNoML3$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML3$condName)])
stripchart(df_ReducAkylNoML3$CysModPerc[df_ReducAkylNoML3$Rep==3]~df_ReducAkylNoML3$condName[df_ReducAkylNoML3$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML3$condName)])
pvalues<-pValueCalSort(df_ReducAkylNoML3$CysModPerc,df_ReducAkylNoML3$condName)
stars<-pvalueStars(pvalues)
arrows(x0=seq(1,4)[!is.na(stars)],x1=seq(2,5)[!is.na(stars)],y0=1.01*sort(Mean_CysModPerc,decreasing = T)[-5][!is.na(stars)],y1=1.01*sort(Mean_CysModPerc,decreasing = T)[-5][!is.na(stars)],code=3,angle=90,length = 0)
text(x=(seq(1,4)[!is.na(stars)]+seq(2,5)[!is.na(stars)])/2,y=1.011*sort(Mean_CysModPerc,decreasing = T)[-5][!is.na(stars)],labels=stars[!is.na(stars)])

stripchart(df_ReducAkylNoML4$CysPeptidePerc[df_ReducAkylNoML4$Rep==1]~df_ReducAkylNoML4$condName[df_ReducAkylNoML4$Rep==1],vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(0,15),col=condCol[levels(df_ReducAkylNoML4$condName)],ylab="CysPeptidePerc")
sem4<-aggregate(df_ReducAkylNoML[,4], list(df_ReducAkylNoML4$condName), FUN=sd)[,2]
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_CysPeptidePerc,decreasing = T)-sem4,
         y1=sort(Mean_CysPeptidePerc,decreasing = T)+sem4)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_CysPeptidePerc,decreasing = T)-sem4,
         y1=sort(Mean_CysPeptidePerc,decreasing = T)-sem4)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_CysPeptidePerc,decreasing = T)+sem4,
         y1=sort(Mean_CysPeptidePerc,decreasing = T)+sem4)
stripchart(df_ReducAkylNoML4$CysPeptidePerc[df_ReducAkylNoML4$Rep==1]~df_ReducAkylNoML4$condName[df_ReducAkylNoML4$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(0,15),col=condCol[levels(df_ReducAkylNoML4$condName)],ylab="CysPeptidePerc")
stripchart(df_ReducAkylNoML4$CysPeptidePerc[df_ReducAkylNoML4$Rep==2]~df_ReducAkylNoML4$condName[df_ReducAkylNoML4$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML4$condName)])
stripchart(df_ReducAkylNoML4$CysPeptidePerc[df_ReducAkylNoML4$Rep==3]~df_ReducAkylNoML4$condName[df_ReducAkylNoML4$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML4$condName)])
segments(x0=seq(1,5)-0.1,x1=seq(1,5)+0.1,y0=sort(Mean_CysPeptidePerc,decreasing = T),y1=sort(Mean_CysPeptidePerc,decreasing = T))
pvalues<-pValueCalSort(df_ReducAkylNoML4$CysPeptidePerc,df_ReducAkylNoML4$condName)
stars<-pvalueStars(pvalues)
arrows(x0=seq(1,4)[!is.na(stars)],x1=seq(2,5)[!is.na(stars)],y0=1.1*sort(Mean_CysPeptidePerc,decreasing = T)[-5][!is.na(stars)],y1=1.1*sort(Mean_CysPeptidePerc,decreasing = T)[-5][!is.na(stars)],code=3,angle=90,length = 0)
text(x=(seq(1,4)[!is.na(stars)]+seq(2,5)[!is.na(stars)])/2,y=1.14*sort(Mean_CysPeptidePerc,decreasing = T)[-5][!is.na(stars)],labels=stars[!is.na(stars)])



stripchart(df_ReducAkylNoML1$NumPeptide[df_ReducAkylNoML1$Rep==1]~df_ReducAkylNoML1$condName[df_ReducAkylNoML1$Rep==1],vertical=T,method = "jitter",las=2,cex=0,pch=15,ylim=c(0,30000),col=condCol[levels(df_ReducAkylNoML1$condName)]
           ,ylab="NumPeptide")
segments(x0=seq(1,5)-0.1,x1=seq(1,5)+0.1,y0=sort(Mean_NumPeptide,decreasing = T),y1=sort(Mean_NumPeptide,decreasing = T))
sem1<-aggregate(df_ReducAkylNoML[,1], list(df_ReducAkylNoML1$condName), FUN=sd)[,2]
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_NumPeptide,decreasing = T)-sem1,
         y1=sort(Mean_NumPeptide,decreasing = T)+sem1)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_NumPeptide,decreasing = T)-sem1,
         y1=sort(Mean_NumPeptide,decreasing = T)-sem1)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_NumPeptide,decreasing = T)+sem1,
         y1=sort(Mean_NumPeptide,decreasing = T)+sem1)
stripchart(df_ReducAkylNoML1$NumPeptide[df_ReducAkylNoML1$Rep==1]~df_ReducAkylNoML1$condName[df_ReducAkylNoML1$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(0,30000),col=condCol[levels(df_ReducAkylNoML1$condName)]
           ,ylab="NumPeptide")
stripchart(df_ReducAkylNoML1$NumPeptide[df_ReducAkylNoML1$Rep==2]~df_ReducAkylNoML1$condName[df_ReducAkylNoML1$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML1$condName)])
stripchart(df_ReducAkylNoML1$NumPeptide[df_ReducAkylNoML1$Rep==3]~df_ReducAkylNoML1$condName[df_ReducAkylNoML1$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML1$condName)])
pvalues<-pValueCalSort(df_ReducAkylNoML1$NumPeptide,df_ReducAkylNoML1$condName)
stars<-pvalueStars(pvalues)
arrows(x0=seq(1,4)[!is.na(stars)],x1=seq(2,5)[!is.na(stars)],y0=1.1*sort(Mean_NumPeptide,decreasing = T)[-5][!is.na(stars)],y1=1.1*sort(Mean_NumPeptide,decreasing = T)[-5][!is.na(stars)],code=3,angle=90,length = 0)
text(x=(seq(1,4)[!is.na(stars)]+seq(2,5)[!is.na(stars)])/2,y=1.14*sort(Mean_NumPeptide,decreasing = T)[-5][!is.na(stars)],labels=stars[!is.na(stars)])


stripchart(df_ReducAkylNoML2$NumProtein[df_ReducAkylNoML2$Rep==1]~df_ReducAkylNoML2$condName[df_ReducAkylNoML2$Rep==1],vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(0,3500),col=condCol[levels(df_ReducAkylNoML2$condName)],ylab="NumProtein")
segments(x0=seq(1,5)-0.1,x1=seq(1,5)+0.1,y0=sort(Mean_NumProtein,decreasing = T),y1=sort(Mean_NumProtein,decreasing = T))
sem2<-aggregate(df_ReducAkylNoML[,2], list(df_ReducAkylNoML2$condName), FUN=sd)[,2]
segments(x0=seq(1,5),x1=seq(1,5),y0=sort(Mean_NumProtein,decreasing = T)-sem2,
         y1=sort(Mean_NumProtein,decreasing = T)+sem2)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_NumProtein,decreasing = T)-sem2,
         y1=sort(Mean_NumProtein,decreasing = T)-sem2)
segments(x0=seq(1,5)-0.05,x1=seq(1,5)+0.05,y0=sort(Mean_NumProtein,decreasing = T)+sem2,
         y1=sort(Mean_NumProtein,decreasing = T)+sem2)
stripchart(df_ReducAkylNoML2$NumProtein[df_ReducAkylNoML2$Rep==1]~df_ReducAkylNoML2$condName[df_ReducAkylNoML2$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(0,3500),col=condCol[levels(df_ReducAkylNoML2$condName)],ylab="NumProtein")
stripchart(df_ReducAkylNoML2$NumProtein[df_ReducAkylNoML2$Rep==2]~df_ReducAkylNoML2$condName[df_ReducAkylNoML2$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML2$condName)])
stripchart(df_ReducAkylNoML2$NumProtein[df_ReducAkylNoML2$Rep==3]~df_ReducAkylNoML2$condName[df_ReducAkylNoML2$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML2$condName)])
pvalues<-pValueCalSort(df_ReducAkylNoML2$NumProtein,df_ReducAkylNoML2$condName)
stars<-pvalueStars(pvalues)
if(sum(is.na(stars))!=length(stars)){
  arrows(x0=seq(1,4)[!is.na(stars)],x1=seq(2,5)[!is.na(stars)],y0=1.1*sort(Mean_NumProtein,decreasing = T)[-5][!is.na(stars)],y1=1.1*sort(Mean_NumProtein,decreasing = T)[-5][!is.na(stars)],code=3,angle=90,length = 0)
  text(x=(seq(1,4)[!is.na(stars)]+seq(2,5)[!is.na(stars)])/2,y=1.14*sort(Mean_NumProtein,decreasing = T)[-5][!is.na(stars)],labels=stars[!is.na(stars)])
}
dev.off()
}

# Homogenization####
preMat_filename = '2_Homogenization_report.pr_matrix.tsv'
proMat_filename = '2_Homogenization_report.pg_matrix.tsv'

preMatIn<-read.table(preMat_filename,stringsAsFactors = F,sep='\t',header = T)
proMatIn<-read.table(proMat_filename,stringsAsFactors = F,sep='\t',header = T)

removeNames<-function(colnamein){
  return(paste(unlist(strsplit(unlist(strsplit(colnamein,"\\."))[13],"_"))[c(2,5)],collapse = "_"))
}
preCond<-sapply(colnames(preMatIn[,11:ncol(preMatIn)]),removeNames)
pepAll<-unique(preMatIn$Stripped.Sequence)
df_PepID<-as.data.frame(matrix(0,nrow=length(pepAll),ncol=length(preCond)))
colnames(df_PepID)<-preCond
row.names(df_PepID)<-pepAll
for(j in 11:ncol(preMatIn)){
  pepAllSel<-unique(preMatIn$Stripped.Sequence[!is.na(preMatIn[,j])])
  df_PepID[pepAllSel,j-10]<-1
}
Homo_pepNum<-apply(df_PepID,2,sum)

proMat<-proMatIn[,6:ncol(proMatIn)]
row.names(proMat)<-proMatIn$Protein.Group
colnames(proMat)<-sapply(colnames(proMat),removeNames)
df_proID<-proMat
df_proID[!is.na(df_proID)]<-1
df_proID[is.na(df_proID)]<-0

Homo_proNum<-apply(df_proID,2,sum)
Homo_Num<-as.data.frame(cbind(Homo_pepNum,Homo_proNum))
Homo_Num$GelName<-factor(sapply(sapply(row.names(Homo_Num),strsplit,"_"),"[[",1))
Homo_col<-hcl.colors(4, palette = "Warm", alpha = 0.9)
names(Homo_col)<-as.vector(levels(Homo_Num$GelName))


Mean_NumPeptide<-aggregate(Homo_Num$Homo_pepNum, list(Homo_Num$GelName), FUN=mean)[,2]
Mean_NumProtein<-aggregate(Homo_Num$Homo_proNum, list(Homo_Num$GelName), FUN=mean)[,2]

{
pdf("2_Homogenization_pvalue.pdf",width=5,height=2.71)
par(mfrow=c(1,2),mar=c(4,4,1,1))
stripchart(Homo_Num$Homo_pepNum~Homo_Num$GelName,col=Homo_col,ylab="NumPeptide",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,40000),cex=0)
segments(x0=seq(1,4)-0.1,x1=seq(1,4)+0.1,y0=Mean_NumPeptide,y1=Mean_NumPeptide)
sem1<-aggregate(Homo_Num[,1], list(Homo_Num$GelName), FUN=sd)[,2]
segments(x0=seq(1,4),x1=seq(1,4),y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide+sem1)
segments(x0=seq(1,4)-0.05,x1=seq(1,4)+0.05,y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide-sem1)
segments(x0=seq(1,4)-0.05,x1=seq(1,4)+0.05,y0=Mean_NumPeptide+sem1,
         y1=Mean_NumPeptide+sem1)
stripchart(Homo_Num$Homo_pepNum~Homo_Num$GelName,col=Homo_col,ylab="NumPeptide",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,40000),add=T,cex=0.7)
pvalues<-pValueCal(Homo_Num$Homo_pepNum,Homo_Num$GelName)
stars<-pvalueStars(pvalues)
combSel<-combn(4,2)[,!is.na(stars)]
arrows(x0=combSel[1,],x1=combSel[2,],y0=1.5*(Mean_NumPeptide[combSel[1,]]+Mean_NumPeptide[combSel[2,]])/2,y1=1.5*(Mean_NumPeptide[combSel[1,]]+Mean_NumPeptide[combSel[2,]])/2,code=3,angle=90,length = 0)
text(x=(combSel[1,]+combSel[2,])/2,y=1.51*(Mean_NumPeptide[combSel[1,]]+Mean_NumPeptide[combSel[2,]])/2,labels=stars[!is.na(stars)])



stripchart(Homo_Num$Homo_proNum~Homo_Num$GelName,col=Homo_col,ylab="NumProtein",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,5000),cex=0)
segments(x0=seq(1,4)-0.1,x1=seq(1,4)+0.1,y0=Mean_NumProtein,y1=Mean_NumProtein)
sem1<-aggregate(Homo_Num[,2], list(Homo_Num$GelName), FUN=sd)[,2]
segments(x0=seq(1,4),x1=seq(1,4),y0=Mean_NumProtein-sem1,
         y1=Mean_NumProtein+sem1)
segments(x0=seq(1,4)-0.05,x1=seq(1,4)+0.05,y0=Mean_NumProtein-sem1,
         y1=Mean_NumProtein-sem1)
segments(x0=seq(1,4)-0.05,x1=seq(1,4)+0.05,y0=Mean_NumProtein+sem1,
         y1=Mean_NumProtein+sem1)
stripchart(Homo_Num$Homo_proNum~Homo_Num$GelName,col=Homo_col,ylab="NumProtein",vertical=T,method='jitter',pch=19
           ,las=2,ylim=c(0,5000),add=T,cex=0.7)
pvalues<-pValueCal(Homo_Num$Homo_proNum,Homo_Num$GelName)
stars<-pvalueStars(pvalues)
combSel<-combn(4,2)[,!is.na(stars)]
arrows(x0=combSel[1,],x1=combSel[2,],y0=1.25*(Mean_NumProtein[combSel[1,]]+Mean_NumProtein[combSel[2,]])/2,y1=1.25*(Mean_NumProtein[combSel[1,]]+Mean_NumProtein[combSel[2,]])/2,code=3,angle=90,length = 0)
text(x=(combSel[1,]+combSel[2,])/2,y=1.26*(Mean_NumProtein[combSel[1,]]+Mean_NumProtein[combSel[2,]])/2,labels=stars[!is.na(stars)])
dev.off()
}
# LC####
preMat_filename = '3_LC_gradients_report.pr_matrix.tsv'
proMat_filename = '3_LC_gradients_report.pg_matrix.tsv'

preMatIn<-read.table(preMat_filename,stringsAsFactors = F,sep='\t',header = T)
proMatIn<-read.table(proMat_filename,stringsAsFactors = F,sep='\t',header = T)

removeNames<-function(colnamein){
  return(paste(unlist(strsplit(unlist(strsplit(colnamein,"\\."))[13],"_"))[c(3,5)],collapse = "_"))
}
preCond<-sapply(colnames(preMatIn[,11:ncol(preMatIn)]),removeNames)
pepAll<-unique(preMatIn$Stripped.Sequence)
df_PepID<-as.data.frame(matrix(0,nrow=length(pepAll),ncol=length(preCond)))
colnames(df_PepID)<-preCond
row.names(df_PepID)<-pepAll
for(j in 11:ncol(preMatIn)){
  pepAllSel<-unique(preMatIn$Stripped.Sequence[!is.na(preMatIn[,j])])
  df_PepID[pepAllSel,j-10]<-1
}
LC_pepNum<-apply(df_PepID,2,sum)

proMat<-proMatIn[,6:ncol(proMatIn)]
row.names(proMat)<-proMatIn$Protein.Group
colnames(proMat)<-sapply(colnames(proMat),removeNames)
df_proID<-proMat
df_proID[!is.na(df_proID)]<-1
df_proID[is.na(df_proID)]<-0

LC_proNum<-apply(df_proID,2,sum)
LC_Num<-as.data.frame(cbind(LC_pepNum,LC_proNum))
LC_Num$GradLength<-factor(sapply(sapply(row.names(LC_Num),strsplit,"_"),"[[",1))
LC_col<-hcl.colors(4, palette = "Cold", alpha = 0.9)
names(LC_col)<-as.vector(levels(LC_Num$GradLength))

LC_Num$MeanPeak<-c(2.46,2.51,2.53,3.47,3.55,3.4)

Mean_NumPeptide<-aggregate(LC_Num$LC_pepNum, list(LC_Num$GradLength), FUN=mean)[,2]
Mean_NumProtein<-aggregate(LC_Num$LC_proNum, list(LC_Num$GradLength), FUN=mean)[,2]
Mean_MeanPeak<-aggregate(LC_Num$MeanPeak, list(LC_Num$GradLength), FUN=mean)[,2]

sem1<-aggregate(LC_Num[,1], list(LC_Num$GradLength), FUN=sd)[,2]
sem2<-aggregate(LC_Num[,2], list(LC_Num$GradLength), FUN=sd)[,2]
sem3<-aggregate(LC_Num[,4], list(LC_Num$GradLength), FUN=sd)[,2]


{
  pdf("3_LC_gradients_pvalue.pdf",width=5,height=2.71)
  par(mfrow=c(1,3),mar=c(4,4,1,1))
  stripchart(LC_Num$LC_pepNum~LC_Num$GradLength,col=LC_col,ylab="NumPeptide",vertical=T,method='jitter',pch=19
             ,las=2,ylim=c(0,50000),cex=0,at=c(1,2),xlim=c(0.5,2.5))
  segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_NumPeptide,y1=Mean_NumPeptide)
  segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_NumPeptide-sem1,
           y1=Mean_NumPeptide+sem1)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumPeptide-sem1,
           y1=Mean_NumPeptide-sem1)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumPeptide+sem1,
           y1=Mean_NumPeptide+sem1)
  stripchart(LC_Num$LC_pepNum~LC_Num$GradLength,col=LC_col,ylab="NumPeptide",vertical=T,method='jitter',pch=19
             ,las=2,ylim=c(0,50000),add=T,cex=0.7,at=c(1,2),xlim=c(0,3))
  pvalues<-pValueCal(LC_Num$LC_pepNum,LC_Num$GradLength)
  stars<-pvalueStars(pvalues)
  arrows(x0=1,x1=2,y0=1.2*(Mean_NumPeptide[1]+Mean_NumPeptide[2])/2,y1=1.2*(Mean_NumPeptide[1]+Mean_NumPeptide[2])/2,code=3,angle=90,length = 0)
  text(x=1.5,y=1.21*(Mean_NumPeptide[1]+Mean_NumPeptide[2])/2,labels=stars[!is.na(stars)])
  
  
  stripchart(LC_Num$LC_proNum~LC_Num$GradLength,col=LC_col,ylab="NumProtein",vertical=T,method='jitter',pch=19
             ,las=2,ylim=c(0,6000),cex=0,at=c(1,2),xlim=c(0.5,2.5))
  segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_NumProtein,y1=Mean_NumProtein)
  segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_NumProtein-sem2,
           y1=Mean_NumProtein+sem2)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumProtein-sem2,
           y1=Mean_NumProtein-sem2)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_NumProtein+sem2,
           y1=Mean_NumProtein+sem2)
  stripchart(LC_Num$LC_proNum~LC_Num$GradLength,col=LC_col,ylab="NumProtein",vertical=T,method='jitter',pch=19
             ,las=2,ylim=c(0,6000),add=T,cex=0.7,at=c(1,2))
  pvalues<-pValueCal(LC_Num$LC_proNum,LC_Num$GradLength)
  stars<-pvalueStars(pvalues)
  arrows(x0=1,x1=2,y0=1.1*(Mean_NumProtein[1]+Mean_NumProtein[2])/2,y1=1.1*(Mean_NumProtein[1]+Mean_NumProtein[2])/2,code=3,angle=90,length = 0)
  text(x=1.5,y=1.11*(Mean_NumProtein[1]+Mean_NumProtein[2])/2,labels=stars[!is.na(stars)])
  
  
  stripchart(LC_Num$MeanPeak~LC_Num$GradLength,col=LC_col,ylab="MeanPeakFWHM",vertical=T,method='jitter',pch=19,las=2,ylim=c(0,4),cex=0,at=c(1,2),xlim=c(0.5,2.5))
  segments(x0=seq(1,2)-0.1,x1=seq(1,2)+0.1,y0=Mean_MeanPeak,y1=Mean_MeanPeak)
  segments(x0=seq(1,2),x1=seq(1,2),y0=Mean_MeanPeak-sem3,
           y1=Mean_MeanPeak+sem3)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MeanPeak-sem3,
           y1=Mean_MeanPeak-sem3)
  segments(x0=seq(1,2)-0.05,x1=seq(1,2)+0.05,y0=Mean_MeanPeak+sem3,
           y1=Mean_MeanPeak+sem3)
  stripchart(LC_Num$MeanPeak~LC_Num$GradLength,col=LC_col,ylab="MeanPeakFWHM",vertical=T,method='jitter',pch=19,las=2,ylim=c(0,4),add=T,cex=0.7,at=c(1,2))
  pvalues<-pValueCal(LC_Num$MeanPeak,LC_Num$GradLength)
  stars<-pvalueStars(pvalues)
  arrows(x0=1,x1=2,y0=1.24*(Mean_MeanPeak[1]+Mean_MeanPeak[2])/2,y1=1.24*(Mean_MeanPeak[1]+Mean_MeanPeak[2])/2,code=3,angle=90,length = 0)
  text(x=1.5,y=1.25*(Mean_MeanPeak[1]+Mean_MeanPeak[2])/2,labels=stars[!is.na(stars)])
  dev.off()
}

##SFigure-Performance of reduction and alkylation at different stages using mouse brain slices.
fileNames<-c("20211226_peptides_red_alk.csv","20211226_proteins_red_alk.csv","20211226_ratios_of_peptide_containing_cysteine.csv","20211226_carbamidomethyl_modification_rate.csv","20211226_missed cleavage rates_red_alk.csv")

# ReductionAlkylations####
df_ReducAkyl<-as.data.frame(matrix(NA,nrow=3*3,ncol=4))
colnames(df_ReducAkyl)<-c("NumPeptide","NumProtein","CysModPerc","CysPeptidePerc")
df_ReducAkyl$condName<-c(rep("Before embedding",3),rep("After homogenisation",3),rep("During in-gel digestion",3))
df_ReducAkyl$Rep<-rep(c(1, 2, 3),3)

#pepitde number
readPep<-read.csv("20211226_peptides_red_alk.csv",stringsAsFactors = F)[,1:3]
readPro<-read.csv("20211226_proteins_red_alk.csv",stringsAsFactors = F,header = T)[,1:3]
readCysModPerc<-read.csv("20211226_carbamidomethyl_modification_rate.csv",stringsAsFactors = F,header = T)[,1:3]
readCysPeptidePerc<-read.csv("20211226_ratios_of_peptide_containing_cysteine.csv",stringsAsFactors = F,header = T)[,1:3]


df_ReducAkyl$NumPeptide<-unlist(readPep)
df_ReducAkyl$NumProtein<-unlist(readPro)
df_ReducAkyl$CysModPerc<-unlist(readCysModPerc)
df_ReducAkyl$CysPeptidePerc<-unlist(readCysPeptidePerc)

df_ReducAkyl$condName<-as.factor(df_ReducAkyl$condName)
par(mfrow=c(1,4),mar=c(8,4,2,2))
stripchart(NumPeptide~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,25000))
stripchart(NumProtein~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,3400))
stripchart(CysModPerc~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(50,100))
stripchart(CysPeptidePerc~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,15))

df_ReducAkylNoML<-df_ReducAkyl
df_ReducAkylNoML$condName<-factor(df_ReducAkylNoML$condName,levels=condLevel)

Mean_NumPeptide<-aggregate(df_ReducAkylNoML$NumPeptide, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_NumProtein<-aggregate(df_ReducAkylNoML$NumProtein, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_CysModPerc<-aggregate(df_ReducAkylNoML$CysModPerc, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_CysPeptidePerc<-aggregate(df_ReducAkylNoML$CysPeptidePerc, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
names(Mean_NumPeptide)<-names(Mean_NumProtein)<-names(Mean_CysModPerc)<-names(Mean_CysPeptidePerc)<-c("Before embedding","After homogenisation","During in-gel digestion")
condLevel<-c("Before embedding","After homogenisation","During in-gel digestion")

condCol<-hcl.colors(3, "Set 2")
names(condCol)<-condLevel
pdf("S_ReductionAlkylations_Gel.pdf",width=8.27,height=2.71)
par(mfrow=c(1,4),mar=c(8,4,1,1))
stripchart(df_ReducAkylNoML$CysModPerc~df_ReducAkylNoML$condName,vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(50,100),col=condCol,ylab="CysModPerc")
segments(x0=seq(1,3)-0.1,x1=seq(1,3)+0.1,y0=Mean_CysModPerc,y1=Mean_CysModPerc)
sem3<-aggregate(df_ReducAkylNoML3$CysModPerc, list(df_ReducAkylNoML3$condName), FUN=sd)[,2]
segments(x0=seq(1,3),x1=seq(1,3),y0=Mean_CysModPerc-sem3,
         y1=Mean_CysModPerc+sem3)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_CysModPerc-sem3,
         y1=Mean_CysModPerc-sem3)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_CysModPerc+sem3,
         y1=Mean_CysModPerc+sem3)
stripchart(df_ReducAkylNoML$CysModPerc[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,add=T,cex=0.7,pch=15,ylim=c(50,100),col=condCol,ylab="CysModPerc")
stripchart(df_ReducAkylNoML$CysModPerc[df_ReducAkylNoML$Rep==2]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol)
stripchart(df_ReducAkylNoML$CysModPerc[df_ReducAkylNoML$Rep==3]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==3],vertical=T,method = "jitter",las=2,add=T,cex=0.7,pch=17,col=condCol)


stripchart(df_ReducAkylNoML$CysPeptidePerc[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(0,8),col=condCol[levels(df_ReducAkylNoML$condName)],ylab="CysPeptidePerc")
sem4<-aggregate(df_ReducAkylNoML[,4], list(df_ReducAkylNoML$condName), FUN=sd)[,2]
segments(x0=seq(1,3),x1=seq(1,3),y0=Mean_CysPeptidePerc-sem4,
         y1=Mean_CysPeptidePerc+sem4)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_CysPeptidePerc-sem4,
         y1=Mean_CysPeptidePerc-sem4)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_CysPeptidePerc+sem4,
         y1=Mean_CysPeptidePerc+sem4)
stripchart(df_ReducAkylNoML$CysPeptidePerc[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(0,5),col=condCol[levels(df_ReducAkylNoML$condName)],ylab="CysPeptidePerc")
stripchart(df_ReducAkylNoML$CysPeptidePerc[df_ReducAkylNoML$Rep==2]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML$condName)])
stripchart(df_ReducAkylNoML$CysPeptidePerc[df_ReducAkylNoML$Rep==3]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML$condName)])
segments(x0=seq(1,3)-0.1,x1=seq(1,3)+0.1,y0=Mean_CysPeptidePerc,y1=Mean_CysPeptidePerc)


stripchart(df_ReducAkylNoML$NumPeptide[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,cex=0,pch=15,ylim=c(0,25000),col=condCol[levels(df_ReducAkylNoML$condName)]
           ,ylab="NumPeptide")
segments(x0=seq(1,3)-0.1,x1=seq(1,3)+0.1,y0=Mean_NumPeptide,y1=Mean_NumPeptide)
sem1<-aggregate(df_ReducAkylNoML[,1], list(df_ReducAkylNoML$condName), FUN=sd)[,2]
segments(x0=seq(1,3),x1=seq(1,3),y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide+sem1)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide-sem1)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_NumPeptide+sem1,
         y1=Mean_NumPeptide+sem1)
stripchart(df_ReducAkylNoML$NumPeptide[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(0,25000),col=condCol[levels(df_ReducAkylNoML$condName)]
           ,ylab="NumPeptide")
stripchart(df_ReducAkylNoML$NumPeptide[df_ReducAkylNoML$Rep==2]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML$condName)])
stripchart(df_ReducAkylNoML$NumPeptide[df_ReducAkylNoML$Rep==3]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML$condName)])


stripchart(df_ReducAkylNoML$NumProtein[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(0,3400),col=condCol[levels(df_ReducAkylNoML$condName)],ylab="NumProtein")
segments(x0=seq(1,3)-0.1,x1=seq(1,3)+0.1,y0=Mean_NumProtein,y1=Mean_NumProtein)
sem2<-aggregate(df_ReducAkylNoML[,2], list(df_ReducAkylNoML$condName), FUN=sd)[,2]
segments(x0=seq(1,3),x1=seq(1,3),y0=Mean_NumProtein-sem2,
         y1=Mean_NumProtein+sem2)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_NumProtein-sem2,
         y1=Mean_NumProtein-sem2)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_NumProtein+sem2,
         y1=Mean_NumProtein+sem2)
stripchart(df_ReducAkylNoML$NumProtein[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=0.7,pch=15,ylim=c(0,3400),col=condCol[levels(df_ReducAkylNoML$condName)],ylab="NumProtein")
stripchart(df_ReducAkylNoML$NumProtein[df_ReducAkylNoML$Rep==2]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol[levels(df_ReducAkylNoML$condName)])
stripchart(df_ReducAkylNoML$NumProtein[df_ReducAkylNoML$Rep==3]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=17,col=condCol[levels(df_ReducAkylNoML$condName)])
dev.off()

#peptide length
readCysPepLen<-read.csv("20211126_peptide length_red_alk.csv",stringsAsFactors = F,header = T,row.names = 1)[,1:4]
pdf("S_Gel_PepLen.pdf",height = 2.7,width=4.1)
plot(readCysPepLen$Peptide.length-0.2,readCysPepLen$Before.embedding,col=condCol[1],type='h',lwd=1.7,xlab="Peptide lengths",ylab="#",axes=F
     ,ylim=c(0,8000))
lines(readCysPepLen$Peptide.length,readCysPepLen$After.homogenisation,col=condCol[2],type='h',lwd=1.7)
lines(readCysPepLen$Peptide.length+0.2,readCysPepLen$During.in.gel.digestion,col=condCol[3],type='h',lwd=1.7)

axis(side=1,at=seq(7,50,3),labels = seq(7,50,3))
axis(side=2,at=seq(0,8000,2000),labels = seq(0,8000,2000),las=2)
legend("topright",legend=c("Red/alk before embedding (whole gel)","Red/alk after homogenisation (whole gel)","Red/alk during in-gel digestion (single punch)"
),lwd=1.7,col=condCol)
dev.off()

pdf("S_Gel_PepLen_curve.pdf",height = 2.7,width=4.1)

plot(readCysPepLen$Peptide.length,readCysPepLen$Before.embedding,col=condCol[1],type='l',lwd=8,xlab="Peptide lengths",ylab="#",axes=F,
     ylim=c(0,8000))
lines(readCysPepLen$Peptide.length,readCysPepLen$After.homogenisation,col=condCol[2],type='l',lwd=8)
lines(readCysPepLen$Peptide.length+0.2,readCysPepLen$During.in.gel.digestion,col=condCol[3],type='l',lwd=8)

axis(side=1,at=seq(7,50,3),labels = seq(7,50,3))
axis(side=2,at=seq(0,8000,2000),labels = seq(0,8000,2000),las=2)
dev.off()


#mis-cleavages 
readMisClv<-read.csv("20211226_missed cleavage rates_red_alk.csv",stringsAsFactors = F,header = T,row.names = 1)[,1:9]
misClvRate<-(readMisClv[2,]+readMisClv[3,])
df_ReducAkylNoML$MisClvRate<-as.numeric(misClvRate)

Mean_MisClvRate<-aggregate(df_ReducAkylNoML$MisClvRate, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
names(Mean_MisClvRate)<-c("Before embedding","After homogenisation","During in-gel digestion")

pdf("S_ReductionAlkylationsMisClv_Gel.pdf",width=2.71,height=3.11)
par(mfrow=c(1,1),mar=c(5,4,1,1))
stripchart(df_ReducAkylNoML$MisClvRate~df_ReducAkylNoML$condName,vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(0,25),col=condCol,ylab="Missed cleavage rate(%)")
segments(x0=seq(1,3)-0.1,x1=seq(1,3)+0.1,y0=Mean_MisClvRate,y1=Mean_MisClvRate)
sem3<-aggregate(df_ReducAkylNoML$MisClvRate, list(df_ReducAkylNoML$condName), FUN=sd)[,2]
segments(x0=seq(1,3),x1=seq(1,3),y0=Mean_MisClvRate-sem3,
         y1=Mean_MisClvRate+sem3)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_MisClvRate-sem3,
         y1=Mean_MisClvRate-sem3)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_MisClvRate+sem3,
         y1=Mean_MisClvRate+sem3)
stripchart(df_ReducAkylNoML$MisClvRate[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,add=T,cex=0.7,pch=15,ylim=c(50,100),col=condCol,ylab="CysModPerc")
stripchart(df_ReducAkylNoML$MisClvRate[df_ReducAkylNoML$Rep==2]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=0.7,pch=16,col=condCol)
stripchart(df_ReducAkylNoML$MisClvRate[df_ReducAkylNoML$Rep==3]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==3],vertical=T,method = "jitter",las=2,add=T,cex=0.7,pch=17,col=condCol)
dev.off()

##SFigure-Performance of different enzyme digestion conditions using mouse brain slices.
# Trypsin####
df_ReducAkyl<-as.data.frame(matrix(NA,nrow=3*3,ncol=2))
colnames(df_ReducAkyl)<-c("NumPeptide","NumProtein")
df_ReducAkyl$condName<-c(rep("Trypsin 4+12h",3),rep("LysC + Trypsin",3),rep("Trypsin 3h",3))
df_ReducAkyl$Rep<-rep(c(1, 2, 3),3)

readPep<-read.csv("20211226_peptides_trp.csv",stringsAsFactors = F)[,1:3]
readPro<-read.csv("20211226_proteins_trp.csv",stringsAsFactors = F,header = T)[,1:3]


df_ReducAkyl$NumPeptide<-unlist(readPep)
df_ReducAkyl$NumProtein<-unlist(readPro)

df_ReducAkyl$condName<-as.factor(df_ReducAkyl$condName)
par(mfrow=c(1,2),mar=c(8,4,2,2))
stripchart(NumPeptide~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,25000))
stripchart(NumProtein~condName,data=df_ReducAkyl,vertical=T,method = "jitter",las=2
           ,ylim=c(0,3400))

df_ReducAkylNoML<-df_ReducAkyl
condLevel<-c("Trypsin 4+12h","LysC + Trypsin","Trypsin 3h")
df_ReducAkylNoML$condName<-factor(df_ReducAkylNoML$condName,levels=condLevel)

Mean_NumPeptide<-aggregate(df_ReducAkylNoML$NumPeptide, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
Mean_NumProtein<-aggregate(df_ReducAkylNoML$NumProtein, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
names(Mean_NumPeptide)<-names(Mean_NumProtein)<-c("Trypsin 4+12h","LysC + Trypsin","Trypsin 3h")

condCol<-hcl.colors(3, "Set 2")
names(condCol)<-condLevel
pdf("S_Trypsin_Gel.pdf",width=8.27,height=5)
par(mfrow=c(1,2),mar=c(8,4,2,2))

stripchart(df_ReducAkylNoML$NumPeptide[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,cex=0,pch=15,ylim=c(0,25000),col=condCol[levels(df_ReducAkylNoML$condName)]
           ,ylab="NumPeptide")
segments(x0=seq(1,3)-0.1,x1=seq(1,3)+0.1,y0=Mean_NumPeptide,y1=Mean_NumPeptide)
sem1<-aggregate(df_ReducAkylNoML[,1], list(df_ReducAkylNoML$condName), FUN=sd)[,2]
segments(x0=seq(1,3),x1=seq(1,3),y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide+sem1)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_NumPeptide-sem1,
         y1=Mean_NumPeptide-sem1)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_NumPeptide+sem1,
         y1=Mean_NumPeptide+sem1)
stripchart(df_ReducAkylNoML$NumPeptide[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=1.7,pch=15,ylim=c(0,25000),col=condCol[levels(df_ReducAkylNoML$condName)]
           ,ylab="NumPeptide")
stripchart(df_ReducAkylNoML$NumPeptide[df_ReducAkylNoML$Rep==2]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=1.7,pch=16,col=condCol[levels(df_ReducAkylNoML$condName)])
stripchart(df_ReducAkylNoML$NumPeptide[df_ReducAkylNoML$Rep==3]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=1.7,pch=17,col=condCol[levels(df_ReducAkylNoML$condName)])


stripchart(df_ReducAkylNoML$NumProtein[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(0,3400),col=condCol[levels(df_ReducAkylNoML$condName)],ylab="NumProtein")
segments(x0=seq(1,3)-0.1,x1=seq(1,3)+0.1,y0=Mean_NumProtein,y1=Mean_NumProtein)
sem2<-aggregate(df_ReducAkylNoML[,2], list(df_ReducAkylNoML$condName), FUN=sd)[,2]
segments(x0=seq(1,3),x1=seq(1,3),y0=Mean_NumProtein-sem2,
         y1=Mean_NumProtein+sem2)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_NumProtein-sem2,
         y1=Mean_NumProtein-sem2)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_NumProtein+sem2,
         y1=Mean_NumProtein+sem2)
stripchart(df_ReducAkylNoML$NumProtein[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,add=T
           ,cex=1.7,pch=15,ylim=c(0,3400),col=condCol[levels(df_ReducAkylNoML$condName)],ylab="NumProtein")
stripchart(df_ReducAkylNoML$NumProtein[df_ReducAkylNoML$Rep==2]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=1.7,pch=16,col=condCol[levels(df_ReducAkylNoML$condName)])
stripchart(df_ReducAkylNoML$NumProtein[df_ReducAkylNoML$Rep==3]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==3],vertical=T,method = "jitter",las=2
           ,add=T,cex=1.7,pch=17,col=condCol[levels(df_ReducAkylNoML$condName)])
dev.off()

#peptide length
readCysPepLen<-read.csv("20211126_peptide length_tryp.csv",stringsAsFactors = F,header = T,row.names = 1)[,1:4]
pdf("S_TryP_curve.pdf",height = 2.7,width=4.1)

plot(readCysPepLen$Peptide.length,readCysPepLen$Trypsin.4...12h,col=condCol[1],type='l',lwd=8,xlab="Peptide lengths",ylab="#",axes=F,
     ylim=c(0,8000))
lines(readCysPepLen$Peptide.length,readCysPepLen$LysC...Trypsin,col=condCol[2],type='l',lwd=8)
lines(readCysPepLen$Peptide.length+0.2,readCysPepLen$Trypsin.3h,col=condCol[3],type='l',lwd=8)

axis(side=1,at=seq(7,50,3),labels = seq(7,50,3))
axis(side=2,at=seq(0,8000,2000),labels = seq(0,8000,2000),las=2)
legend("topright", legend = c("Trypsin 4+12h", "LysC + Trypsin", "Trypsin 3h"), lwd = 1.7, col = condCol)
text("topright", pos = c(1,1), "Enzyme Treatment", cex = 0.8, adj = c(1, 0))
dev.off()

#mis-cleavages 
readMisClv<-read.csv("20211226_missed cleavage rates_trp.csv",stringsAsFactors = F,header = T,row.names = 1)[,1:9]
misClvRate<-(readMisClv[2,]+readMisClv[3,])
df_ReducAkylNoML$MisClvRate<-as.numeric(misClvRate)

Mean_MisClvRate<-aggregate(df_ReducAkylNoML$MisClvRate, list(df_ReducAkylNoML$condName), FUN=mean)[,2]
names(Mean_MisClvRate)<-c("Trypsin 4+12h","LysC + Trypsin","Trypsin 3h")

pdf("S_TrypsinMisClv_Gel.pdf",width=2.71,height=3.11)
par(mfrow=c(1,1),mar=c(5,4,1,1))
stripchart(df_ReducAkylNoML$MisClvRate~df_ReducAkylNoML$condName,vertical=T,method = "jitter",las=2
           ,cex=0,pch=15,ylim=c(0,25),col=condCol,ylab="Missed cleavage rate(%)")
segments(x0=seq(1,3)-0.1,x1=seq(1,3)+0.1,y0=Mean_MisClvRate,y1=Mean_MisClvRate)
sem3<-aggregate(df_ReducAkylNoML$MisClvRate, list(df_ReducAkylNoML$condName), FUN=sd)[,2]
segments(x0=seq(1,3),x1=seq(1,3),y0=Mean_MisClvRate-sem3,
         y1=Mean_MisClvRate+sem3)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_MisClvRate-sem3,
         y1=Mean_MisClvRate-sem3)
segments(x0=seq(1,3)-0.05,x1=seq(1,3)+0.05,y0=Mean_MisClvRate+sem3,
         y1=Mean_MisClvRate+sem3)
stripchart(df_ReducAkylNoML$MisClvRate[df_ReducAkylNoML$Rep==1]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==1],vertical=T,method = "jitter",las=2,add=T,cex=1.7,pch=15,ylim=c(50,100),col=condCol,ylab="CysModPerc")
stripchart(df_ReducAkylNoML$MisClvRate[df_ReducAkylNoML$Rep==2]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==2],vertical=T,method = "jitter",las=2
           ,add=T,cex=1.7,pch=16,col=condCol)
stripchart(df_ReducAkylNoML$MisClvRate[df_ReducAkylNoML$Rep==3]~df_ReducAkylNoML$condName[df_ReducAkylNoML$Rep==3],vertical=T,method = "jitter",las=2,add=T,cex=1.7,pch=17,col=condCol)
dev.off()