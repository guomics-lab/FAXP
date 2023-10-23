# TimeComp<-read.csv("Figure5/Figure5B.csv",stringsAsFactors = F)
# The data is located in the Figure 5 folder under the "All_matrix" folder.
TimeCompOd<-read.csv("Figure5/Figure5B_od.csv",stringsAsFactors = F)

TimeCompV1<-TimeCompOd[TimeCompOd$V==1,]
TimeCompV2<-TimeCompOd[TimeCompOd$V==2,]

TimeV1<-TimeCompV1[,3]-TimeCompV1[,2]
names(TimeV1)<-TimeCompV1$X
TimeV2<-TimeCompV2[,3]-TimeCompV2[,2]
names(TimeV2)<-TimeCompV2$X

par(mfrow=c(3,1),mar=c(4,20,1,1))
barplot(as.matrix(TimeV1),horiz = T,xlim=c(0,60))
barplot(as.matrix(TimeV2),horiz = T,xlim=c(0,60))

# barplot(as.matrix(TimeV1),horiz = T,xlim=c(0,70),beside = T)
# barplot(as.matrix(TimeV2),horiz = T,xlim=c(0,70),beside = T)
TimeV<-TimeCompOd[,3]-TimeCompOd[,2]
names(TimeV)<-paste(TimeCompOd$X,TimeCompOd$V,sep="_v")
barplot(TimeV[length(TimeV):1],horiz = T,xlim=c(0,60),las=2)


colAll<-hcl.colors(nrow(df_task), "Dynamic")
TimeV<-TimeCompOd[,3]-TimeCompOd[,2]
names(TimeV)<-paste(TimeCompOd$X,TimeCompOd$V,sep="_v")

task<-unique(TimeCompOd$X)
df_task<-matrix(NA,nrow=length(task),ncol=2)
row.names(df_task)<-task
colnames(df_task)<-c("v1","v2")
df_task[names(TimeV1),1]<-TimeV1
df_task[names(TimeV2),2]<-TimeV2

TimeV1P<-df_task[,1]
TimeV2P<-df_task[,2]
TimeV1P[is.na(TimeV1P)]<-0
TimeV2P[is.na(TimeV2P)]<-0


pdf("F5B_TimeComp_FZ20230329.pdf",width=12,height = 8)
par(mfrow=c(1,2))
par(mar=c(5,4,12,5))
barplot(as.matrix(TimeV1P[nrow(df_task):1]),horiz = T,xlim=c(0,60),col=colAll,space=0.01,border = 0.1)
par(mar=c(5,5,12,4))
barplot(as.matrix(TimeV2P[nrow(df_task):1]),horiz = T,xlim=c(0,60),col=colAll,space=0.01,border = 0.1)
dev.off()

par(mar=c(2,4,1,5))
barplot((df_task[nrow(df_task):1,1]),horiz = T,xlim=c(60,0),beside = T,las=2,col=colAll,space=0.01,border = 0.1)
par(mar=c(2,5,1,4))
barplot((df_task[nrow(df_task):1,2]),horiz = T,xlim=c(0,60),beside = T,las=2,col=colAll,space=0.01,border = 0.1)
dev.off()


pdf("F5B_TimeCompSegV1_FZ20230329.pdf",width=12,height = 12)
par(mfcol=c(nrow(df_task)+1,1),mar=c(0,1,0,3))
for(i in 1:nrow(df_task)){
  colAllSel<-rep(0,nrow(df_task))
  colAllSel[i]<-colAll[i]
  barplot(as.matrix(TimeV1P[nrow(df_task):1]),horiz = T,xlim=c(0,60),col=colAllSel,space=0.01,border = 0.1,axes=F)
}
par(mar=c(4,1,0,3))
barplot(as.matrix(TimeV1P[nrow(df_task):1]),horiz = T,xlim=c(0,60),col=colAllSel,space=0.01,border = 0.1,axes=T)
dev.off()

pdf("F5B_TimeCompSegV2_FZ20230329.pdf",width=12,height = 12)
par(mfcol=c(nrow(df_task)+1,1),mar=c(0,1,0,3))
for(i in 1:nrow(df_task)){
  colAllSel<-rep(0,nrow(df_task))
  colAllSel[i]<-colAll[i]
  barplot(as.matrix(TimeV2P[nrow(df_task):1]),horiz = T,xlim=c(0,60),col=colAllSel,space=0.01,border = 0.1,axes=F)
}
par(mar=c(4,1,0,3))
barplot(as.matrix(TimeV2P[nrow(df_task):1]),horiz = T,xlim=c(0,60),col=colAllSel,space=0.01,border = 0.1,axes=T)

dev.off()