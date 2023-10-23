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


proMatNoPool<-read.csv("proMatNoPool_FZ20230322.csv",stringsAsFactors = F,row.names = 1)
proMat<-proMatNoPool[!grepl(";",row.names(proMatNoPool)),]
colnames(proMat)<-batchDesign2[colnames(proMatNoPool),"NewID"]
quantile(apply(is.na(proMat),2,sum)/nrow(proMat))
quantile(apply(is.na(proMat),1,sum)/ncol(proMat))
quantile(apply(is.na(proMat),1,sum)/ncol(proMat),probs = seq(0, 1, 0.1))

#Inter variational analysis####
CoV<-function(vectorIN){#vectorIN=c(1,2,3,NA)
  return(sd(vectorIN,na.rm=T)/mean(vectorIN,na.rm=T))
  # return(var(vectorIN,na.rm=T)/mean(vectorIN,na.rm=T))
}
#interpatient
proMat_N<-proMat[,grepl("N",colnames(proMat))]
proMat_L<-proMat[,grepl("L",colnames(proMat))]
proMat_H<-proMat[,grepl("H",colnames(proMat))]
proMat_C<-proMat[,grepl("C",colnames(proMat))]

InterPat_all<-apply(proMat,1,CoV)
InterPat_N<-apply(proMat_N,1,CoV)
InterPat_L<-apply(proMat_L,1,CoV)
InterPat_H<-apply(proMat_H,1,CoV)
InterPat_C<-apply(proMat_C,1,CoV)

InterPatList<-list(InterPat_all,InterPat_N,InterPat_L,InterPat_H,InterPat_C)
names(InterPatList)<-c("all","N","L","H","C")

proMat_NA<-na.omit(proMat)
proMat_N_NA<-na.omit(proMat_N)
proMat_L_NA<-na.omit(proMat_L)
proMat_H_NA<-na.omit(proMat_H)
proMat_C_NA<-na.omit(proMat_C)

InterPat_all_NA<-apply(proMat_NA,1,CoV)
InterPat_N_NA<-apply(proMat_N_NA,1,CoV)
InterPat_L_NA<-apply(proMat_L_NA,1,CoV)
InterPat_H_NA<-apply(proMat_H_NA,1,CoV)
InterPat_C_NA<-apply(proMat_C_NA,1,CoV)

InterPatListNA<-list(InterPat_all_NA,InterPat_N_NA,InterPat_L_NA,InterPat_H_NA,InterPat_C_NA)
names(InterPatListNA)<-c("all","N","L","H","C")

par(mfrow=c(1,2))
vioplot::vioplot(InterPatList,col=0,ylim=c(0,6),main="InterPatient CoV\nwith missing")
boxplot(InterPatList,add=T,outline=F)
vioplot::vioplot(InterPatListNA,col=0,ylim=c(0,6),main="InterPatient CoV\nno missing")
boxplot(InterPatListNA,add=T,outline=F)

#interslide
proMat_P1<-proMat[,grepl("P1",colnames(proMat))]
proMat_P1_N<-proMat[,grepl("N",colnames(proMat))&grepl("P1",colnames(proMat))]
proMat_P1_L<-proMat[,grepl("L",colnames(proMat))&grepl("P1",colnames(proMat))]
proMat_P1_H<-proMat[,grepl("H",colnames(proMat))&grepl("P1",colnames(proMat))]
proMat_P1_PC<-proMat[,grepl("PC",colnames(proMat))&grepl("P1",colnames(proMat))]
proMat_P1_CC<-proMat[,grepl("CC",colnames(proMat))&grepl("P1",colnames(proMat))]

InterSlide_P1<-apply(proMat_P1,1,CoV)
InterSlide_N<-apply(proMat_P1_N,1,CoV)
InterSlide_L<-apply(proMat_P1_L,1,CoV)
InterSlide_H<-apply(proMat_P1_H,1,CoV)
InterSlide_PC<-apply(proMat_P1_PC,1,CoV)
InterSlide_CC<-apply(proMat_P1_CC,1,CoV)

InterSlideList<-list(InterSlide_P1,InterSlide_N,InterSlide_L,
                   InterSlide_H,InterSlide_PC,InterSlide_CC)
names(InterSlideList)<-c("P1","N","L","H","PC","CC")

par(mfrow=c(1,2))
vioplot::vioplot(InterPatList,col=0,ylim=c(0,6),main="InterPatient CoV")
boxplot(InterPatList,add=T,outline=F)
vioplot::vioplot(InterSlideList,col=0,ylim=c(0,6),main="InterSlide CoV")
boxplot(InterSlideList,add=T,outline=F)

#Intratissue
proMat_P1S1<-proMat[,grepl("P1S1",colnames(proMat))]
proMat_P1S1_N<-proMat[,grepl("N",colnames(proMat))&grepl("P1S1",colnames(proMat))]
proMat_P1S1_L<-proMat[,grepl("L",colnames(proMat))&grepl("P1S1",colnames(proMat))]
proMat_P1S1_H<-proMat[,grepl("H",colnames(proMat))&grepl("P1S1",colnames(proMat))]
proMat_P1S1_PC<-proMat[,grepl("PC",colnames(proMat))&grepl("P1S1",colnames(proMat))]
proMat_P1S1_CC<-proMat[,grepl("CC",colnames(proMat))&grepl("P1S1",colnames(proMat))]

proMat_P1S2<-proMat[,grepl("P1S2",colnames(proMat))]
proMat_P1S2_N<-proMat[,grepl("N",colnames(proMat))&grepl("P1S2",colnames(proMat))]
proMat_P1S2_L<-proMat[,grepl("L",colnames(proMat))&grepl("P1S2",colnames(proMat))]
proMat_P1S2_H<-proMat[,grepl("H",colnames(proMat))&grepl("P1S2",colnames(proMat))]
proMat_P1S2_PC<-proMat[,grepl("PC",colnames(proMat))&grepl("P1S2",colnames(proMat))]
proMat_P1S2_CC<-proMat[,grepl("CC",colnames(proMat))&grepl("P1S2",colnames(proMat))]

proMat_P1S3<-proMat[,grepl("P1S3",colnames(proMat))]
proMat_P1S3_N<-proMat[,grepl("N",colnames(proMat))&grepl("P1S3",colnames(proMat))]
proMat_P1S3_L<-proMat[,grepl("L",colnames(proMat))&grepl("P1S3",colnames(proMat))]
proMat_P1S3_H<-proMat[,grepl("H",colnames(proMat))&grepl("P1S3",colnames(proMat))]
proMat_P1S3_PC<-proMat[,grepl("PC",colnames(proMat))&grepl("P1S3",colnames(proMat))]
proMat_P1S3_CC<-proMat[,grepl("CC",colnames(proMat))&grepl("P1S3",colnames(proMat))]

proMat_P2S1<-proMat[,grepl("P2S1",colnames(proMat))]
proMat_P2S1_L<-proMat[,grepl("L",colnames(proMat))&grepl("P2S1",colnames(proMat))]
proMat_P2S1_H<-proMat[,grepl("H",colnames(proMat))&grepl("P2S1",colnames(proMat))]
proMat_P2S1_C<-proMat[,grepl("C",colnames(proMat))&grepl("P2S1",colnames(proMat))]

proMat_P3S1<-proMat[,grepl("P3S1",colnames(proMat))]
proMat_P3S1_N<-proMat[,grepl("N",colnames(proMat))&grepl("P3S1",colnames(proMat))]
proMat_P3S1_L<-proMat[,grepl("L",colnames(proMat))&grepl("P3S1",colnames(proMat))]
proMat_P3S1_C<-proMat[,grepl("C",colnames(proMat))&grepl("P3S1",colnames(proMat))]

IntraTissue_P1S1<-apply(proMat_P1S1,1,CoV)
IntraTissue_P1S1_N<-apply(proMat_P1S1_N,1,CoV)
IntraTissue_P1S1_L<-apply(proMat_P1S1_L,1,CoV)
IntraTissue_P1S1_H<-apply(proMat_P1S1_H,1,CoV)
IntraTissue_P1S1_PC<-apply(proMat_P1S1_PC,1,CoV)
IntraTissue_P1S1_CC<-apply(proMat_P1S1_CC,1,CoV)
IntraTissue_P1S1List<-list(IntraTissue_P1S1,IntraTissue_P1S1_N,IntraTissue_P1S1_L,
                           IntraTissue_P1S1_H,IntraTissue_P1S1_PC,IntraTissue_P1S1_CC)
names(IntraTissue_P1S1List)<-c("P1S1","N","L","H","PC","CC")

IntraTissue_P1S2<-apply(proMat_P1S2,1,CoV)
IntraTissue_P1S2_N<-apply(proMat_P1S2_N,1,CoV)
IntraTissue_P1S2_L<-apply(proMat_P1S2_L,1,CoV)
IntraTissue_P1S2_H<-apply(proMat_P1S2_H,1,CoV)
IntraTissue_P1S2_PC<-apply(proMat_P1S2_PC,1,CoV)
IntraTissue_P1S2_CC<-apply(proMat_P1S2_CC,1,CoV)
IntraTissue_P1S2List<-list(IntraTissue_P1S2,IntraTissue_P1S2_N,IntraTissue_P1S2_L,
                           IntraTissue_P1S2_H,IntraTissue_P1S2_PC,IntraTissue_P1S2_CC)
names(IntraTissue_P1S2List)<-c("P1S2","N","L","H","PC","CC")

IntraTissue_P1S3<-apply(proMat_P1S3,1,CoV)
IntraTissue_P1S3_N<-apply(proMat_P1S3_N,1,CoV)
IntraTissue_P1S3_L<-apply(proMat_P1S3_L,1,CoV)
IntraTissue_P1S3_H<-apply(proMat_P1S3_H,1,CoV)
IntraTissue_P1S3_PC<-apply(proMat_P1S3_PC,1,CoV)
IntraTissue_P1S3_CC<-apply(proMat_P1S3_CC,1,CoV)
IntraTissue_P1S3List<-list(IntraTissue_P1S3,IntraTissue_P1S3_N,IntraTissue_P1S3_L,
                           IntraTissue_P1S3_H,IntraTissue_P1S3_PC,IntraTissue_P1S3_CC)
names(IntraTissue_P1S3List)<-c("P1S3","N","L","H","PC","CC")

IntraTissue_P2S1<-apply(proMat_P2S1,1,CoV)
IntraTissue_P2S1_L<-apply(proMat_P2S1_L,1,CoV)
IntraTissue_P2S1_H<-apply(proMat_P2S1_H,1,CoV)
IntraTissue_P2S1_C<-apply(proMat_P2S1_C,1,CoV)
IntraTissue_P2S1List<-list(IntraTissue_P2S1,IntraTissue_P2S1_L,IntraTissue_P2S1_H,IntraTissue_P2S1_C)
names(IntraTissue_P2S1List)<-c("P2S1","L","H","C")

IntraTissue_P3S1<-apply(proMat_P3S1,1,CoV)
IntraTissue_P3S1_N<-apply(proMat_P3S1_N,1,CoV)
IntraTissue_P3S1_L<-apply(proMat_P3S1_L,1,CoV)
IntraTissue_P3S1_C<-apply(proMat_P3S1_C,1,CoV)
IntraTissue_P3S1List<-list(IntraTissue_P3S1,IntraTissue_P3S1_N,IntraTissue_P3S1_L,IntraTissue_P3S1_C)
names(IntraTissue_P3S1List)<-c("P3S1","N","L","C")

pdf("InterIntraCoV.pdf",width=21,height=3)
par(mfrow=c(1,7))
vioplot::vioplot(InterPatList,col=2,ylim=c(0,1),main="InterPatient CoV")
boxplot(InterPatList,add=T,outline=F)
vioplot::vioplot(InterSlideList,col=3,ylim=c(0,1),main="InterSlide CoV")
boxplot(InterSlideList,add=T,outline=F)
vioplot::vioplot(IntraTissue_P1S1List,col=4,ylim=c(0,1),main="IntraTissue P1S1 CoV")
boxplot(IntraTissue_P1S1List,add=T,outline=F)
vioplot::vioplot(IntraTissue_P1S2List,col=5,ylim=c(0,1),main="IntraTissue P1S2 CoV")
boxplot(IntraTissue_P1S2List,add=T,outline=F)
vioplot::vioplot(IntraTissue_P1S3List,col=6,ylim=c(0,1),main="IntraTissue P1S3 CoV")
boxplot(IntraTissue_P1S3List,add=T,outline=F)
vioplot::vioplot(IntraTissue_P2S1List,col=7,ylim=c(0,1),main="IntraTissue P2S1 CoV")
boxplot(IntraTissue_P2S1List,add=T,outline=F)
vioplot::vioplot(IntraTissue_P3S1List,col=8,ylim=c(0,1),main="IntraTissue P3S1 CoV")
boxplot(IntraTissue_P3S1List,add=T,outline=F)
dev.off()


pdf("InterIntraCoV_groupedbox.pdf",width=12,height=4)
inspace=8
stInd<-seq(0,by=1,length.out=7)
par(mfrow=c(1,1))
boxplot(InterPatList,outline=F,at=seq(stInd[1],by=inspace,length.out=5),col=2,xlim=c(0,55))
boxplot(InterSlideList,outline=F,at=seq(stInd[2],by=inspace,length.out=6),col=3,add=T)
boxplot(IntraTissue_P1S1List,outline=F,at=seq(stInd[3],by=inspace,length.out=6),col=4,add=T)
boxplot(IntraTissue_P1S2List,outline=F,at=seq(stInd[4],by=inspace,length.out=6),col=5,add=T)
boxplot(IntraTissue_P1S3List,outline=F,at=seq(stInd[5],by=inspace,length.out=6),col=6,add=T)
boxplot(IntraTissue_P2S1List,outline=F,at=seq(stInd[6],by=inspace,length.out=6)[c(1,3,4,5)],col=7,add=T)
boxplot(IntraTissue_P3S1List,outline=F,at=seq(stInd[7],by=inspace,length.out=6)[c(1,2,3,5)],col=8,add=T)
legend("topright",fill=seq(2,8),legend = c("InterPat","InterSlide","IntraTissueP1S1","IntraTissueP1S2","IntraTissueP1S3"
                                           ,"IntraTissueP2S1","IntraTissueP3S1"))
dev.off()

#no na
IntraTissue_P1S1_NA<-apply(na.omit(proMat_P1S1),1,CoV)
IntraTissue_P1S1_N_NA<-apply(na.omit(proMat_P1S1_N),1,CoV)
IntraTissue_P1S1_L_NA<-apply(na.omit(proMat_P1S1_L),1,CoV)
IntraTissue_P1S1_H_NA<-apply(na.omit(proMat_P1S1_H),1,CoV)
IntraTissue_P1S1_PC_NA<-apply(na.omit(proMat_P1S1_PC),1,CoV)
IntraTissue_P1S1_CC_NA<-apply(na.omit(proMat_P1S1_CC),1,CoV)
IntraTissue_P1S1List_NA<-list(IntraTissue_P1S1_NA,IntraTissue_P1S1_N_NA,IntraTissue_P1S1_L_NA,
                           IntraTissue_P1S1_H_NA,IntraTissue_P1S1_PC_NA,IntraTissue_P1S1_CC_NA)
names(IntraTissue_P1S1List_NA)<-c("P1S1","N","L","H","PC","CC")

IntraTissue_P1S2_NA<-apply(na.omit(proMat_P1S2),1,CoV)
IntraTissue_P1S2_N_NA<-apply(na.omit(proMat_P1S2_N),1,CoV)
IntraTissue_P1S2_L_NA<-apply(na.omit(proMat_P1S2_L),1,CoV)
IntraTissue_P1S2_H_NA<-apply(na.omit(proMat_P1S2_H),1,CoV)
IntraTissue_P1S2_PC_NA<-apply(na.omit(proMat_P1S2_PC),1,CoV)
IntraTissue_P1S2_CC_NA<-apply(na.omit(proMat_P1S2_CC),1,CoV)
IntraTissue_P1S2List_NA<-list(IntraTissue_P1S2_NA,IntraTissue_P1S2_N_NA,IntraTissue_P1S2_L_NA,
                           IntraTissue_P1S2_H_NA,IntraTissue_P1S2_PC_NA,IntraTissue_P1S2_CC_NA)
names(IntraTissue_P1S2List_NA)<-c("P1S2","N","L","H","PC","CC")

IntraTissue_P1S3_NA<-apply(na.omit(proMat_P1S3),1,CoV)
IntraTissue_P1S3_N_NA<-apply(na.omit(proMat_P1S3_N),1,CoV)
IntraTissue_P1S3_L_NA<-apply(na.omit(proMat_P1S3_L),1,CoV)
IntraTissue_P1S3_H_NA<-apply(na.omit(proMat_P1S3_H),1,CoV)
IntraTissue_P1S3_PC_NA<-apply(na.omit(proMat_P1S3_PC),1,CoV)
IntraTissue_P1S3_CC_NA<-apply(na.omit(proMat_P1S3_CC),1,CoV)
IntraTissue_P1S3List_NA<-list(IntraTissue_P1S3_NA,IntraTissue_P1S3_N_NA,IntraTissue_P1S3_L_NA,
                           IntraTissue_P1S3_H_NA,IntraTissue_P1S3_PC_NA,IntraTissue_P1S3_CC_NA)
names(IntraTissue_P1S3List_NA)<-c("P1S3","N","L","H","PC","CC")

par(mfrow=c(1,5))
vioplot::vioplot(InterPatList,col=0,ylim=c(0,1),main="InterPatient CoV")
boxplot(InterPatList,add=T,outline=F)
vioplot::vioplot(InterSlideList,col=0,ylim=c(0,1),main="InterSlide CoV")
boxplot(InterSlideList,add=T,outline=F)
vioplot::vioplot(IntraTissue_P1S1List_NA,col=0,ylim=c(0,1),main="IntraTissue P1S1 CoV")
boxplot(IntraTissue_P1S1List_NA,add=T,outline=F)
vioplot::vioplot(IntraTissue_P1S2List_NA,col=0,ylim=c(0,1),main="IntraTissue P1S2 CoV")
boxplot(IntraTissue_P1S2List_NA,add=T,outline=F)
vioplot::vioplot(IntraTissue_P1S3List_NA,col=0,ylim=c(0,1),main="IntraTissue P1S3 CoV")
boxplot(IntraTissue_P1S3List_NA,add=T,outline=F)

inspace=6
stInd<-seq(0,by=1,length.out=5)
boxplot(InterPatList,outline=F,at=seq(stInd[1],by=inspace,length.out=5),col=2,xlim=c(0,40))
boxplot(InterSlideList,outline=F,at=seq(stInd[2],by=inspace,length.out=6),col=3,add=T)
boxplot(IntraTissue_P1S1List_NA,outline=F,at=seq(stInd[3],by=inspace,length.out=6),col=4,add=T)
boxplot(IntraTissue_P1S2List_NA,outline=F,at=seq(stInd[4],by=inspace,length.out=6),col=5,add=T)
boxplot(IntraTissue_P1S3List_NA,outline=F,at=seq(stInd[5],by=inspace,length.out=6),col=6,add=T)
legend("topright",fill=seq(2,6),legend = c("InterPat","InterSlide","IntraTissueS1","IntraTissueS2","IntraTissueS3"))

#InterPunch
proMat_P1_CC1<-proMat[,grepl("P1",colnames(proMat))&grepl("CC1",colnames(proMat))]
proMat_P1_CC2<-proMat[,grepl("P1",colnames(proMat))&grepl("CC2",colnames(proMat))]
proMat_P1_CC3<-proMat[,grepl("P1",colnames(proMat))&grepl("CC3",colnames(proMat))]
proMat_P1_CC4<-proMat[,grepl("P1",colnames(proMat))&grepl("CC4",colnames(proMat))]
proMat_P1_CC5<-proMat[,grepl("P1",colnames(proMat))&grepl("CC5",colnames(proMat))]

IntraPunch_P1_CC1<-apply(proMat_P1_CC1,1,CoV)
IntraPunch_P1_CC2<-apply(proMat_P1_CC2,1,CoV)
IntraPunch_P1_CC3<-apply(proMat_P1_CC3,1,CoV)
IntraPunch_P1_CC4<-apply(proMat_P1_CC4,1,CoV)
IntraPunch_P1_CC5<-apply(proMat_P1_CC5,1,CoV)
P1_CCList<-list(IntraTissue_P1S1_CC,IntraTissue_P1S2_CC,IntraTissue_P1S3_CC
     ,IntraPunch_P1_CC1,IntraPunch_P1_CC2,IntraPunch_P1_CC3
     ,IntraPunch_P1_CC4,IntraPunch_P1_CC5)
names(P1_CCList)<-c("P1S1_CC","P1S2_CC","P1S3_CC"
                    ,"P1_CC1","P1_CC2","P1_CC3","P1_CC4","P1_CC5")
vioplot::vioplot(P1_CCList,ylim=c(0,1),main="CC IntraTissue vs InterPunch CoV"
                 ,col=0)
boxplot(P1_CCList,add=T,outline=F)


proMat_P1_PC1<-proMat[,grepl("P1",colnames(proMat))&grepl("PC1",colnames(proMat))]
proMat_P1_PC2<-proMat[,grepl("P1",colnames(proMat))&grepl("PC2",colnames(proMat))]
proMat_P1_PC3<-proMat[,grepl("P1",colnames(proMat))&grepl("PC3",colnames(proMat))]
proMat_P1_PC4<-proMat[,grepl("P1",colnames(proMat))&grepl("PC4",colnames(proMat))]
proMat_P1_PC5<-proMat[,grepl("P1",colnames(proMat))&grepl("PC5",colnames(proMat))]
proMat_P1_PC6<-proMat[,grepl("P1",colnames(proMat))&grepl("PC6",colnames(proMat))]

IntraPunch_P1_PC1<-apply(proMat_P1_PC1,1,CoV)
IntraPunch_P1_PC2<-apply(proMat_P1_PC2,1,CoV)
IntraPunch_P1_PC3<-apply(proMat_P1_PC3,1,CoV)
IntraPunch_P1_PC4<-apply(proMat_P1_PC4,1,CoV)
IntraPunch_P1_PC5<-apply(proMat_P1_PC5,1,CoV)
IntraPunch_P1_PC6<-apply(proMat_P1_PC6,1,CoV)

P1_PCList<-list(IntraTissue_P1S1_PC,IntraTissue_P1S2_PC,IntraTissue_P1S3_PC
                ,IntraPunch_P1_PC1,IntraPunch_P1_PC2,IntraPunch_P1_PC3
                ,IntraPunch_P1_PC4,IntraPunch_P1_PC5,IntraPunch_P1_PC6)
names(P1_PCList)<-c("P1S1_PC","P1S2_PC","P1S3_PC"
                    ,"P1_PC1","P1_PC2","P1_PC3","P1_PC4","P1_PC5","P1_PC6")
vioplot::vioplot(P1_PCList,ylim=c(0,1),main="PC IntraTissue vs InterPunch CoV"
                 ,col=0)
boxplot(P1_PCList,add=T,outline=F)

# The following P1_H sections were added by dongzhen on 20230802
proMat_P1_H1<-proMat[,grepl("P1",colnames(proMat))&grepl("H1",colnames(proMat))]
proMat_P1_H2<-proMat[,grepl("P1",colnames(proMat))&grepl("H2",colnames(proMat))]
proMat_P1_H3<-proMat[,grepl("P1",colnames(proMat))&grepl("H3",colnames(proMat))]
proMat_P1_H4<-proMat[,grepl("P1",colnames(proMat))&grepl("H4",colnames(proMat))]
proMat_P1_H5<-proMat[,grepl("P1",colnames(proMat))&grepl("H5",colnames(proMat))]
proMat_P1_H6<-proMat[,grepl("P1",colnames(proMat))&grepl("H6",colnames(proMat))]

IntraPunch_P1_H1<-apply(proMat_P1_H1,1,CoV)
IntraPunch_P1_H2<-apply(proMat_P1_H2,1,CoV)
IntraPunch_P1_H3<-apply(proMat_P1_H3,1,CoV)
IntraPunch_P1_H4<-apply(proMat_P1_H4,1,CoV)
IntraPunch_P1_H5<-apply(proMat_P1_H5,1,CoV)
IntraPunch_P1_H6<-apply(proMat_P1_H6,1,CoV)

P1_HList<-list(IntraTissue_P1S1_H,IntraTissue_P1S2_H,IntraTissue_P1S3_H
                ,IntraPunch_P1_H1,IntraPunch_P1_H2,IntraPunch_P1_H3
                ,IntraPunch_P1_H4,IntraPunch_P1_H5,IntraPunch_P1_H6)
names(P1_HList)<-c("P1S1_H","P1S2_H","P1S3_H"
                    ,"P1_H1","P1_H2","P1_H3","P1_H4","P1_H5","P1_H6")
vioplot::vioplot(P1_HList,ylim=c(0,1),main="H IntraTissue vs InterPunch CoV"
                 ,col=0)
boxplot(P1_HList,add=T,outline=F)


proMat_P1_L1<-proMat[,grepl("P1",colnames(proMat))&grepl("L1",colnames(proMat))]
proMat_P1_L2<-proMat[,grepl("P1",colnames(proMat))&grepl("L2",colnames(proMat))]
proMat_P1_L3<-proMat[,grepl("P1",colnames(proMat))&grepl("L3",colnames(proMat))]
proMat_P1_L4<-proMat[,grepl("P1",colnames(proMat))&grepl("L4",colnames(proMat))]
proMat_P1_L5<-proMat[,grepl("P1",colnames(proMat))&grepl("L5",colnames(proMat))]

IntraPunch_P1_L1<-apply(proMat_P1_L1,1,CoV)
IntraPunch_P1_L2<-apply(proMat_P1_L2,1,CoV)
IntraPunch_P1_L3<-apply(proMat_P1_L3,1,CoV)
IntraPunch_P1_L4<-apply(proMat_P1_L4,1,CoV)
IntraPunch_P1_L5<-apply(proMat_P1_L5,1,CoV)
P1_LList<-list(IntraTissue_P1S1_L,IntraTissue_P1S2_L,IntraTissue_P1S3_L
                ,IntraPunch_P1_L1,IntraPunch_P1_L2,IntraPunch_P1_L3
                ,IntraPunch_P1_L4,IntraPunch_P1_L5)
names(P1_LList)<-c("P1S1_L","P1S2_L","P1S3_L"
                    ,"P1_L1","P1_L2","P1_L3","P1_L4","P1_L5")
vioplot::vioplot(P1_LList,ylim=c(0,1),main="L IntraTissue vs InterPunch CoV"
                 ,col=0)
boxplot(P1_LList,add=T,outline=F)


proMat_P1_N1<-proMat[,grepl("P1",colnames(proMat))&grepl("N1",colnames(proMat))]
proMat_P1_N2<-proMat[,grepl("P1",colnames(proMat))&grepl("N2",colnames(proMat))]
proMat_P1_N3<-proMat[,grepl("P1",colnames(proMat))&grepl("N3",colnames(proMat))]
proMat_P1_N4<-proMat[,grepl("P1",colnames(proMat))&grepl("N4",colnames(proMat))]
proMat_P1_N5<-proMat[,grepl("P1",colnames(proMat))&grepl("N5",colnames(proMat))]
proMat_P1_N6<-proMat[,grepl("P1",colnames(proMat))&grepl("N6",colnames(proMat))]

IntraPunch_P1_N1<-apply(proMat_P1_N1,1,CoV)
IntraPunch_P1_N2<-apply(proMat_P1_N2,1,CoV)
IntraPunch_P1_N3<-apply(proMat_P1_N3,1,CoV)
IntraPunch_P1_N4<-apply(proMat_P1_N4,1,CoV)
IntraPunch_P1_N5<-apply(proMat_P1_N5,1,CoV)
IntraPunch_P1_N6<-apply(proMat_P1_N6,1,CoV)

P1_NList<-list(IntraTissue_P1S1_N,IntraTissue_P1S2_N,IntraTissue_P1S3_N
                ,IntraPunch_P1_N1,IntraPunch_P1_N2,IntraPunch_P1_N3
                ,IntraPunch_P1_N4,IntraPunch_P1_N5,IntraPunch_P1_N6)
names(P1_NList)<-c("P1S1_N","P1S2_N","P1S3_N"
                    ,"P1_N1","P1_N2","P1_N3","P1_N4","P1_N5","P1_N6")
vioplot::vioplot(P1_NList,ylim=c(0,1),main="N IntraTissue vs InterPunch CoV"
                 ,col=0)
boxplot(P1_NList,add=T,outline=F)

ColSubtype<-hcl.colors(6, "viridis", rev = TRUE)[1:5]

pdf("InterTissueVsPunch.pdf",width=20,height=5)
par(mfrow=c(1,5))
vioplot::vioplot(P1_NList,ylim=c(0,1),main="N IntraTissue vs InterPunch CoV"
                 ,col=ColSubtype[1],las=2)
boxplot(P1_NList,add=T,outline=F,axes=F,col=c(rep("white",3),rep(8,length(P1_NList)-3)))
vioplot::vioplot(P1_LList,ylim=c(0,1),main="L IntraTissue vs InterPunch CoV"
                 ,col=ColSubtype[2],las=2)
boxplot(P1_LList,add=T,outline=F,axes=F,col=c(rep("white",3),rep(8,length(P1_LList)-3)))
vioplot::vioplot(P1_HList,ylim=c(0,1),main="H IntraTissue vs InterPunch CoV"
                 ,col=ColSubtype[3],las=2)
boxplot(P1_HList,add=T,outline=F,axes=F,col=c(rep("white",3),rep(8,length(P1_HList)-3)))
vioplot::vioplot(P1_PCList,ylim=c(0,1),main="PC IntraTissue vs InterPunch CoV"
                 ,col=ColSubtype[4],las=2)
boxplot(P1_PCList,add=T,outline=F,axes=F,col=c(rep("white",3),rep(8,length(P1_PCList)-3)))
vioplot::vioplot(P1_CCList,ylim=c(0,1),main="CC IntraTissue vs InterPunch CoV"
                 ,col=ColSubtype[5],las=2)
boxplot(P1_CCList,add=T,outline=F,axes=F,col=c(rep("white",3),rep(8,length(P1_CCList)-3)))
dev.off()

P1_N_ITIP<-list(unlist(P1_NList[1:3]),unlist(P1_NList[4:length(P1_NList)]))
P1_L_ITIP<-list(unlist(P1_LList[1:3]),unlist(P1_LList[4:length(P1_LList)]))
P1_H_ITIP<-list(unlist(P1_HList[1:3]),unlist(P1_HList[4:length(P1_HList)]))
P1_PC_ITIP<-list(unlist(P1_PCList[1:3]),unlist(P1_PCList[4:length(P1_PCList)]))
P1_CC_ITIP<-list(unlist(P1_CCList[1:3]),unlist(P1_CCList[4:length(P1_CCList)]))

pdf("InterTissueVsPunch_2.pdf",width=20,height=5)
par(mfrow=c(1,5))
vioplot(P1_N_ITIP,col=ColSubtype[1],ylim=c(0,1))
bx<-boxplot(P1_N_ITIP,add=T,outline=F,axes=F,col=c("white",'grey'))
text(c(1,2),bx$stats[3,]+0.05,round(bx$stats[3,]*100,2),col=2)
vioplot(P1_L_ITIP,col=ColSubtype[2],ylim=c(0,1))
bx<-boxplot(P1_L_ITIP,add=T,outline=F,axes=F,col=c("white",'grey'))
text(c(1,2),bx$stats[3,]+0.05,round(bx$stats[3,]*100,2),col=2)

vioplot(P1_H_ITIP,col=ColSubtype[3],ylim=c(0,1))
bx<-boxplot(P1_H_ITIP,add=T,outline=F,axes=F,col=c("white",'grey'))
text(c(1,2),bx$stats[3,]+0.05,round(bx$stats[3,]*100,2),col=2)

vioplot(P1_PC_ITIP,col=ColSubtype[4],ylim=c(0,1))
bx<-boxplot(P1_PC_ITIP,add=T,outline=F,axes=F,col=c("white",'grey'))
text(c(1,2),bx$stats[3,]+0.05,round(bx$stats[3,]*100,2),col=2)

vioplot(P1_CC_ITIP,col=ColSubtype[5],ylim=c(0,1))
bx<-boxplot(P1_CC_ITIP,add=T,outline=F,axes=F,col=c("white",'grey'))
text(c(1,2),bx$stats[3,]+0.05,round(bx$stats[3,]*100,2),col=2)
dev.off()

#UMAP####
library(umap)


pchIDs = seq(10,15)
colIDs = hcl.colors(6, "viridis", rev = TRUE)
names(pchIDs)<-names(colIDs)<-c('N','L','H',"C","PC","CC")

NA_thrs=c(0.4,0.5,0.6)
NA_fill_ratios=c(0.7,0.8,0.9)
# pdf("UMAP_NA_thrs0.4_0.6fill0.7_0.9.pdf",width=42,height = 18)
par(mfrow=c(3,7))
for(NA_thr in NA_thrs){
  for(NA_fill_ratio in NA_fill_ratios){
    proMatList<-list(proMat,proMat_P1,proMat_P2S1,proMat_P3S1,proMat_P1S1,proMat_P1S2,proMat_P1S3)
    names(proMatList)<-c("All","P1","P2","P3","P1S1","P1S2","P1S3")
    for (i in 1:length(proMatList)){
      proSel=proMatList[[i]]
      proSel_NA<-proSel[(apply(is.na(proSel),1,sum)/ncol(proSel))<=NA_thr,]
      proSel_NAFill<-proSel_NA
      proSel_NAFill[is.na(proSel_NA)]<-min(proSel_NA,na.rm=T)*NA_fill_ratio
      um<-umap(t(proSel_NAFill))
      subTypesSel<-sapply(colnames(proSel_NAFill),toPlotPar)[3,]
      plot(um$layout[,1],um$layout[,2],xlab='umap1',ylab='umap2'
           ,pch=pchIDs[subTypesSel]
           ,col=colIDs[subTypesSel]
           ,main = paste0(names(proMatList)[i],"NAthr-",NA_thr,"NAfill-",NA_fill_ratio),cex=0.6)
      text(um$layout[,1],um$layout[,2],col=colIDs[subTypesSel],colnames(proSel_NAFill),cex=0.6)
    }
  }
}
# dev.off()

#umap 3D####
library(plot3D)
proMatList<-list(proMat,proMat_P1,proMat_P2S1,proMat_P3S1,proMat_P1S1,proMat_P1S2,proMat_P1S3)
NA_fill_ratio=NA_thr=0.7
names(proMatList)<-c("All","P1","P2","P3","P1S1","P1S2","P1S3")
pdf("Umap2DP1_0.65_0.65.pdf")
for (i in 2:2){#length(proMatList)
  proSel=proMatList[[i]]
  proSel_NA<-proSel[(apply(is.na(proSel),1,sum)/ncol(proSel))<=NA_thr,]
  proSel_NAFill<-proSel_NA
  proSel_NAFill[is.na(proSel_NA)]<-min(proSel_NA,na.rm=T)*NA_fill_ratio
  # um<-umap(t(proSel_NAFill), n_components = 3)
  um<-umap(t(proSel_NAFill), n_components = 2)
  
  subTypesSel<-sapply(colnames(proSel_NAFill),toPlotPar)[3,]
  subTypesSel[subTypesSel=="PC"]<-"C"
  subTypesSel[subTypesSel=="CC"]<-"C"
  
  # plot3d(um$layout[,1], um$layout[,2], um$layout[,3],size=15
            # ,pch=pchIDs[subTypesSel],col=colIDs[subTypesSel]
         # ,xlab='umap1',ylab='umap2',zlab='umap3',aspect=10)
  # 
  plot(um$layout[,1],um$layout[,2],xlab='umap1',ylab='umap2'
       ,pch=19#pchIDs[subTypesSel]
       ,col=colIDs[subTypesSel]
       ,main = paste0(names(proMatList)[i],"NAthr-",NA_thr,"NAfill-",NA_fill_ratio),cex=1.5)
  
  Median_LocY<-Median_LocX<-rep(NA,5)
  names(Median_LocY)<-names(Median_LocX)<-c('N','L','H',"C")
  for(i in 1:length(Median_LocY)){
    Median_LocX[i]<-median(um$layout[,1][grepl(names(Median_LocX)[i],subTypesSel)])
    Median_LocY[i]<-median(um$layout[,2][grepl(names(Median_LocY)[i],subTypesSel)])
  }
  for(i in 1:(length(Median_LocY)-1)){
    arrows(x0=Median_LocX[i],x1=Median_LocX[i+1],
           y0=Median_LocY[i],y1=Median_LocY[i+1],col=colIDs[i],lwd=2)
  }
  # text(um$layout[,1],um$layout[,2],col=colIDs[subTypesSel],colnames(proSel_NAFill),cex=0.6)
}
dev.off()


toPlotPar<-function(NewIDname){
  sss<-unlist(strsplit(NewIDname,""))
  patID<-paste0(sss[1:2],collapse = "")
  textID<-paste0(sss[5:length(sss)],collapse = "")
  subtypeID<-paste0(sss[5:(length(sss)-1)],collapse = "")
  return(c(patID,textID,subtypeID))
}

#hclust####
proMatList<-list(proMat,proMat_P1,proMat_P2S1,proMat_P3S1,proMat_P1S1,proMat_P1S2,proMat_P1S3)
names(proMatList)<-c("All","P1","P2","P3","P1S1","P1S2","P1S3")
NA_thr=0.7
NA_fill_ratio=0.8
pdf("Dendrogram.pdf",width=20,height = 20)
par(mfrow=c(4,1),mar=c(0,0,3,0))#par(mfrow=c(3,1),mar=c(0,1,0,1))
proSel=proMatList[[1]]
proSel_NA<-proSel[(apply(is.na(proSel),1,sum)/ncol(proSel))<=NA_thr,]
proSel_NAFill<-proSel_NA
proSel_NAFill[is.na(proSel_NA)]<-min(proSel_NA,na.rm=T)*NA_fill_ratio
hc <- hclust(dist(t(proSel_NAFill)), "ward.D2")
plot(hc)
orderNames<-colnames(proSel_NAFill)[hc$order]

barPar<-function(NewIDname){
  plotPars<-toPlotPar(NewIDname)
  subType<-plotPars[3]
  PatID<-plotPars[1]
  SlideID<-paste(unlist(strsplit(NewIDname,""))[3:4],collapse = "")
  return(c(PatID,SlideID,subType))
}

GroupNum<-sapply(orderNames,barPar)
GroupNumFactor<-GroupNum
GroupNumFactor[1,]<-as.numeric(as.factor(GroupNum[1,]))*10+100
GroupNumFactor[2,]<-as.numeric(as.factor(GroupNum[2,]))*10+200
GroupNumFactor[3,]<-as.numeric(factor(GroupNum[3,],levels=c('N','L','H',"C","PC","CC")))*10+300
GroupNumMat<-matrix(0,nrow=nrow(GroupNumFactor),ncol=ncol(GroupNumFactor))
GroupNumMat[,]<-as.numeric(GroupNumFactor[,])
par(mar=c(0,2,0,2))
image(t(t(GroupNumMat[3,])),col = hcl.colors(6, "viridis", rev = TRUE),axes=F)
image(t(t(GroupNumMat[2,])),col = hcl.colors(3, "Peach", rev = TRUE),axes=F)
image(t(t(GroupNumMat[1,])),col = hcl.colors(3, "Pastel 1", rev = TRUE),axes=F)
dev.off()

#MultiML####
library(nnet)
library(randomForest)
library(pROC)
set.seed(1717)
NA_thr=0.4;NA_fill_ratio=0.7
subTypes<-c('N','L','H',"C","C","C")
names(subTypes)<-c('N','L','H',"C","PC","CC")
proSel_NA<-proSel[(apply(is.na(proMat),1,sum)/ncol(proSel))<=NA_thr,]
proSel_NAFill<-proSel_NA
proSel_NAFill[is.na(proSel_NA)]<-min(proSel_NA,na.rm=T)*NA_fill_ratio
ml_Mat<-t(proSel_NAFill)
ml_label<-as.factor(as.vector(subTypes[sapply(colnames(proSel_NAFill),barPar)[3,]]))

SampleInd<-sample(nrow(ml_Mat))
ml_Mat_train = ml_Mat[SampleInd[1:(length(SampleInd)*0.7)],]
ml_Mat_test = ml_Mat[SampleInd[(1+length(SampleInd)*0.7):length(SampleInd)],]
ml_label_train = ml_label[SampleInd[1:(length(SampleInd)*0.7)]]
ml_label_test = ml_label[SampleInd[(1+length(SampleInd)*0.7):length(SampleInd)]]

# acc_best=0
for(i in 1:100){
  rf = randomForest(x=ml_Mat_train,y=ml_label_train, ntree = 1717)
  predictions <- as.numeric(predict(rf, ml_Mat_test, type = 'response'))
  predictions_prob <- predict(rf, ml_Mat_test, type = 'prob')
  # roc.multi <- multiclass.roc(ml_label_test, predictions)
  roc.multi_prob <- multiclass.roc(ml_label_test, predictions_prob)
  tabOut<-table(ml_label_test,predict(rf, ml_Mat_test))
  auc = roc.multi_prob$auc
  acc = (tabOut[1,1]+tabOut[2,2]+tabOut[3,3]+tabOut[4,4])/length(ml_label_test)
  if(acc>acc_best){
    rfout<-rf
    print(paste(i,"AUC",round(auc,4),"ACC",round(acc,4)))
    acc_best=acc
    auc_best=auc
  }
  print(i)
}

roc.multi_prob_out <- multiclass.roc(ml_label_test, predict(rfout, ml_Mat_test, type = 'prob'))
tabOut_out<-table(ml_label_test,predict(rfout, ml_Mat_test))

save(SampleInd,rfout,file="rfout_FZ20230323.RDS")

ImpPros<-names(rf$importance[order(rf$importance,decreasing = T),1])[1:100]
write.csv(ImpPros,"ImpPros100.csv",row.names = F)

write.csv(cbind(names(rf$importance[order(rf$importance,decreasing = T),1])[1:100],rf$importance[order(rf$importance,decreasing = T),1][1:100]),"ImpPros100_rfscore.csv",row.names = F)

GnMatch<-read.delim("UniProtGnMatch.tsv",stringsAsFactors = F,sep='\t',row.names = 1,header = T)
for(i in 1:nrow(GnMatch)){
  GnMatch$Gene.Names.Main[i]<-unlist(strsplit(GnMatch$Gene.Names[i]," "))[1]
}
# BiocManager::install("clusterProfiler")
library(pheatmap)
library(RColorBrewer)
toPlotHmap<-proSel_NAFill[ImpPros[1:50],]
row.names(toPlotHmap)<-GnMatch[row.names(toPlotHmap),3]
factorNames<-factor(subTypes[sapply(colnames(toPlotHmap),barPar)[3,]])
ann_colors = list(
  Subtype = c(C = "#7570B3", H = "#E7298A", L = "#66A61E",N="Yellow")
)
annotation_col = data.frame(
  Subtype = factorNames
)
row.names(annotation_col)<-colnames(toPlotHmap)
library(RColorBrewer)
# pdf("Heatmap_CHLN_GN.pdf",width=16,height = 8)
p1<-pheatmap::pheatmap(toPlotHmap[,order(factorNames)],scale='row',cluster_rows=T,cluster_cols=F
         ,color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
         ,annotation_col = annotation_col, annotation_colors = ann_colors[1])
# dev.off()

#spatial images####
install.packages("magick")
library(magick)
HE <- image_read("Figure6_HE_2.png")
HE_scale<-image_scale(HE, "x400")
# img %>% image_convolve('DoG:0,0,2') %>% image_negate()

# image_draw(HE_scale)
img_neg<-image_negate(HE_scale)
image_draw(img_neg)
abline(h = 100, col = 'blue', lwd = '1', lty = "dotted")
abline(v = 300, col = 'blue', lwd = '1', lty = "dotted")


#network####
sub_grp <- cutree(p1$tree_row, k = 4)


mat1<-cor(t(toPlotHmap[names(sub_grp)[sub_grp==1],]))
mat2<-cor(t(toPlotHmap[names(sub_grp)[sub_grp==2],]))
mat3<-cor(t(toPlotHmap[names(sub_grp)[sub_grp==3],]))
mat4<-cor(t(toPlotHmap[names(sub_grp)[sub_grp==4],]))

mat1[mat1<0.7] <- 0
mat2[mat2<0.5] <- 0
mat3[mat3<0.5] <- 0
mat4[mat4<0.5] <- 0

network1 <- graph_from_adjacency_matrix(mat1, weighted=T, mode="undirected", diag=F)
network2 <- graph_from_adjacency_matrix(mat2, weighted=T, mode="undirected", diag=F)
network3 <- graph_from_adjacency_matrix(mat3, weighted=T, mode="undirected", diag=F)
network4 <- graph_from_adjacency_matrix(mat4, weighted=T, mode="undirected", diag=F)

pdf("F6_Network.pdf",width=20,height = 8)
par(mfrow=c(1,4))
plot(network1,node.color="grey80")
plot(network2)
plot(network3)
plot(network4)
dev.off()

#network mfuzz####
sub_grp <- clRanksSel
mat1<-cor(t(proSel_NAFill[names(sub_grp)[sub_grp==0],]))
mat2<-cor(t(proSel_NAFill[names(sub_grp)[sub_grp==1],]))
mat3<-cor(t(proSel_NAFill[names(sub_grp)[sub_grp==2],]))
mat4<-cor(t(proSel_NAFill[names(sub_grp)[sub_grp==3],]))

row.names(mat1)<-colnames(mat1)<-GnMatch[colnames(mat1),3]
row.names(mat2)<-colnames(mat2)<-GnMatch[colnames(mat2),3]
row.names(mat3)<-colnames(mat3)<-GnMatch[colnames(mat3),3]
row.names(mat4)<-colnames(mat4)<-GnMatch[colnames(mat4),3]

mat1[mat1<0.7] <- 0
mat2[mat2<0.5] <- 0
mat3[mat3<0.5] <- 0
mat4[mat4<0.7] <- 0

library(igraph)
network1 <- graph_from_adjacency_matrix(mat1, weighted=T, mode="undirected", diag=F)
network2 <- graph_from_adjacency_matrix(mat2, weighted=T, mode="undirected", diag=F)
network3 <- graph_from_adjacency_matrix(mat3, weighted=T, mode="undirected", diag=F)
network4 <- graph_from_adjacency_matrix(mat4, weighted=T, mode="undirected", diag=F)

pdf("F6_Network_Col10.pdf",width=2,height = 8)
par(mfrow=c(4,1),mar=c(0,0,0,0))
plot(network1, vertex.color = "grey", vertex.label.color="black",edge.color=rainbow(10)
     ,edge.width=edge_attr(network1)$weight*6,layout = layout.circle,vertex.border.width=0.3)
plot(network2, vertex.color = "grey", vertex.label.color="black",edge.color=rainbow(10)
     ,edge.width=edge_attr(network2)$weight*6,layout = layout.circle,vertex.border.width=0.3)
plot(network3, vertex.color = "grey", vertex.label.color="black",edge.color=rainbow(30)
     ,edge.width=edge_attr(network3)$weight*6,layout = layout.circle,vertex.border.width=0.3)
plot(network4, vertex.color = "grey", vertex.label.color="black",edge.color=rainbow(30)
     ,edge.width=edge_attr(network4)$weight*6,layout = layout.circle,vertex.border.width=0.3)
dev.off()
# fviz_cluster(list(data = df, cluster = sub_grp))
#spatial expression####
which(row.names(proMat)=="Q9BQ51")
which(row.names(proMat)=="Q9NZQ7")
which(row.names(proMat)=="Q9BQ51")

sort(colnames(proMat)[grepl("P1S1",colnames(proMat))])
sort(colnames(proMat)[grepl("P1S2",colnames(proMat))])
sort(colnames(proMat)[grepl("P1S3",colnames(proMat))])

SlideLoc<-read.csv("Figure6/SlideLocation.csv",stringsAsFactors = F,header = F,row.names = 1)
SlideLoc_P1S1<-SlideLoc[grepl("P1S1",row.names(SlideLoc)),]
SlideLoc_P1S2<-SlideLoc[grepl("P1S2",row.names(SlideLoc)),]
SlideLoc_P1S3<-SlideLoc[grepl("P1S3",row.names(SlideLoc)),]
# SlideLoc_P1S1_Loc<-SlideLoc_P1S1[1:(nrow(SlideLoc_P1S1)-2),]
# SlideLoc_P1S1_Start<-SlideLoc_P1S1[(nrow(SlideLoc_P1S1)-1),]
# SlideLoc_P1S1_End<-SlideLoc_P1S1[nrow(SlideLoc_P1S1),]
# 
# SlideLoc_P1S1_Loc$x = SlideLoc_P1S1_Loc$V2-SlideLoc_P1S1_Start$V2
# SlideLoc_P1S1_Loc$y = SlideLoc_P1S1_Start$V3-SlideLoc_P1S1_Loc$V3
# plot(SlideLoc_P1S1_Loc$x,SlideLoc_P1S1_Loc$y)
par(mfrow=c(3,1))
plot(SlideLoc_P1S1$V2,SlideLoc_P1S1$V3,ylim=c(200,50))
plot(SlideLoc_P1S2$V2,SlideLoc_P1S2$V3,ylim=c(500,300))
plot(SlideLoc_P1S3$V2,SlideLoc_P1S3$V3,ylim=c(800,600))



colGen<-function(vectorIn,ColIn){
  maxValue<-max(abs(vectorIn),na.rm=T)
  diffValue=maxValue*2/(length(ColIn)-1)
  intIndex = as.integer((vectorIn+maxValue)/diffValue)+1
  return(ColIn[intIndex])
}

# hcl.colors(1000, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE)
hclPatl<-hcl.colors(101, "Blue-Red 3", rev = TRUE)

pdf(paste0("Figure6/","Top9",".pdf"),width=8,height = 12)
for(i in 1:9){
  par(mfrow=c(3,1))
  proSel<-proMat_P1S1[ImpPros[i],]
  proSelZ<-scale(as.numeric(proSel))[,1]
  names(proSelZ)<-colnames(proMat_P1S1)
  plot(SlideLoc_P1S1$V2,SlideLoc_P1S1$V3,ylim=c(200,50),main=ImpPros[i]
       ,col=c(colGen(proSelZ[row.names(SlideLoc_P1S1)[1:28]],hclPatl),0,0),pch=19,cex=3)
  
  proSel<-proMat_P1S2[ImpPros[i],]
  proSelZ<-scale(as.numeric(proSel))[,1]
  names(proSelZ)<-colnames(proMat_P1S2)
  plot(SlideLoc_P1S2$V2,SlideLoc_P1S2$V3,ylim=c(500,300)
       ,col=c(colGen(proSelZ[row.names(SlideLoc_P1S2)[1:28]],hclPatl),0,0),pch=19,cex=3)
  
  proSel<-proMat_P1S3[ImpPros[i],]
  proSelZ<-scale(as.numeric(proSel))[,1]
  names(proSelZ)<-colnames(proMat_P1S3)
  plot(SlideLoc_P1S3$V2,SlideLoc_P1S3$V3,ylim=c(800,600)
       ,col=c(colGen(proSelZ[row.names(SlideLoc_P1S3)[1:28]],hclPatl),0,0),pch=19,cex=3)
}
dev.off()

#three chosen proteins
UniProtSel<-'P01116'
which(row.names(proMat)==UniProtSel)

pdf(paste0("Figure6/","KRAS",".pdf"),width=8,height = 12)
par(mfrow=c(3,1))
proSel<-proMat_P1S1[UniProtSel,]
proSelZ<-scale(as.numeric(proSel))[,1]
names(proSelZ)<-colnames(proMat_P1S1)
plot(SlideLoc_P1S1$V2,SlideLoc_P1S1$V3,ylim=c(200,50),main=UniProtSel
     ,col=c(colGenMax(proSelZ[row.names(SlideLoc_P1S1)[1:28]],hclPatl),0,0),pch=19,cex=3)

proSel<-proMat_P1S2[UniProtSel,]
proSelZ<-scale(as.numeric(proSel))[,1]
names(proSelZ)<-colnames(proMat_P1S2)
plot(SlideLoc_P1S2$V2,SlideLoc_P1S2$V3,ylim=c(500,300)
     ,col=c(colGenMax(proSelZ[row.names(SlideLoc_P1S2)[1:28]],hclPatl),0,0),pch=19,cex=3)

proSel<-proMat_P1S3[UniProtSel,]
proSelZ<-scale(as.numeric(proSel))[,1]
names(proSelZ)<-colnames(proMat_P1S3)
plot(SlideLoc_P1S3$V2,SlideLoc_P1S3$V3,ylim=c(800,600)
     ,col=c(colGenMax(proSelZ[row.names(SlideLoc_P1S3)[1:28]],hclPatl),0,0),pch=19,cex=3)
dev.off()

#CV large cancer
NA_thr=0.1
proMat_P1_N_NA<-proMat_P1_N[(apply(is.na(proMat_P1),1,sum)/ncol(proMat_P1))<=NA_thr,]
proMat_P1_C_NA<-proMat_P1_C[(apply(is.na(proMat_P1),1,sum)/ncol(proMat_P1))<=NA_thr,]
proMat_P1_CC_NA<-proMat_P1_CC[(apply(is.na(proMat_P1),1,sum)/ncol(proMat_P1))<=NA_thr,]
proMat_P1_PC_NA<-proMat_P1_PC[(apply(is.na(proMat_P1),1,sum)/ncol(proMat_P1))<=NA_thr,]

CV_N<-apply(proMat_P1_N_NA,1,sd,na.rm=T)/apply(proMat_P1_N_NA,1,mean,na.rm=T)
CV_C<-apply(proMat_P1_C_NA,1,sd,na.rm=T)/apply(proMat_P1_C_NA,1,mean,na.rm=T)
CV_CC<-apply(proMat_P1_CC_NA,1,sd,na.rm=T)/apply(proMat_P1_CC_NA,1,mean,na.rm=T)
CV_PC<-apply(proMat_P1_PC_NA,1,sd,na.rm=T)/apply(proMat_P1_PC_NA,1,mean,na.rm=T)


UniProtSels<-names(sort(CV_C/CV_N,decreasing = T)[1:9])
UniProtSels<-names(sort((CV_CC+CV_PC)/CV_N,decreasing = T)[1:9])

rank=1
for(UniProtSel in UniProtSels){
  pdf(paste0("Figure6/CV_CCPC_",rank,"_",UniProtSel,".pdf"),width=8,height = 12)
  par(mfrow=c(3,1))
  proSel<-proMat_P1S1[UniProtSel,]
  proSelZ<-scale(as.numeric(proSel))[,1]
  names(proSelZ)<-colnames(proMat_P1S1)
  plot(SlideLoc_P1S1$V2,SlideLoc_P1S1$V3,ylim=c(200,50),main=UniProtSel
       ,col=c(colGenMax(proSelZ[row.names(SlideLoc_P1S1)[1:28]],hclPatl),0,0),pch=19,cex=3)
  
  proSel<-proMat_P1S2[UniProtSel,]
  proSelZ<-scale(as.numeric(proSel))[,1]
  names(proSelZ)<-colnames(proMat_P1S2)
  plot(SlideLoc_P1S2$V2,SlideLoc_P1S2$V3,ylim=c(500,300)
       ,col=c(colGenMax(proSelZ[row.names(SlideLoc_P1S2)[1:28]],hclPatl),0,0),pch=19,cex=3)
  
  proSel<-proMat_P1S3[UniProtSel,]
  proSelZ<-scale(as.numeric(proSel))[,1]
  names(proSelZ)<-colnames(proMat_P1S3)
  plot(SlideLoc_P1S3$V2,SlideLoc_P1S3$V3,ylim=c(800,600)
       ,col=c(colGenMax(proSelZ[row.names(SlideLoc_P1S3)[1:28]],hclPatl),0,0),pch=19,cex=3)
  dev.off()
  rank=rank+1
}


#mfuzz####
library(Mfuzz)
ml_Mat
ml_label
assayMat<-t(ml_Mat)
colnames(assayMat)<-as.vector(ml_label)

assayMatGroup<-matrix(NA,nrow=nrow(assayMat),ncol=4)
colnames(assayMatGroup)<-c('N','L','H',"C")
row.names(assayMatGroup)<-row.names(assayMat)

for(j in 1:ncol(assayMatGroup)){
  assayMatSel<-assayMat[,colnames(assayMat)==colnames(assayMatGroup)[j]]
  assayMatGroup[,j]<-apply(assayMatSel,1,mean)
}
assayMatGroupZ<-assayMatGroup
for(i in 1:nrow(assayMatGroup)){
  assayMatGroupZ[i,]<-scale(assayMatGroup[i,])
}

eset = ExpressionSet(assayData=assayMatGroupZ)

# Soft clustering and visualisation
cln<-4
cl <- mfuzz(eset,c=cln,m=1.15)
pdf("mfuzz_all_4.pdf",width=12,height = 3)
mfuzz.plot(eset,cl=cl,mfrow=c(1,cln),min.mem=0.5,colo=hcl.colors(12, "Blue-Red", rev = TRUE)
           ,time.labels=c('N','L','H',"C"),new.window=F)
dev.off()
#RF features
assayMatGroupZSel<-assayMatGroupZ[ImpPros,]
esetSel = ExpressionSet(assayData=assayMatGroupZSel)
pdf("mfuzz_t100_2.pdf",width=12,height = 3)
cln<-2
cl2 <- mfuzz(esetSel,c=cln,m=1.15)
mfuzz.plot(esetSel,cl=cl2,mfrow=c(1,cln),min.mem=0.5,colo=hcl.colors(12, "Green-Orange", rev = TRUE)
           ,time.labels=c('N','L','H',"C"),new.window=F)
dev.off()

pdf("mfuzz_t100_4.pdf",width=12,height = 3)
cln<-4
cl3 <- mfuzz(esetSel,c=cln,m=1.15)
mfuzz.plot(esetSel,cl=cl3,mfrow=c(1,cln),min.mem=0.5,colo=hcl.colors(12, "Green-Orange", rev = TRUE)
           ,time.labels=c('N','L','H',"C"),new.window=F)
dev.off()

#heatmap mfuzz clusters####
clRanks<-cl3$cluster
clRanks<-clRanks[ImpPros]

clRanks[clRanks==4]=0
clRanksSel<-c(clRanks[clRanks==0][1:10],clRanks[clRanks==1][1:10],
  clRanks[clRanks==2][1:10],clRanks[clRanks==3][1:10])

write.csv(clRanks,"clRanks.csv")

library(pheatmap)
library(RColorBrewer)
toPlotHmap<-proSel_NAFill[names(sort(clRanksSel)),]
row.names(toPlotHmap)<-GnMatch[row.names(toPlotHmap),3]
factorNames<-factor(subTypes[sapply(colnames(toPlotHmap),barPar)[3,]])
ann_colors = list(
  Subtype = c(C = "#008298", H = "#00B28A", L = "#7ED357",N="#FDE333")
)
annotation_col = data.frame(
  Subtype = factorNames
)
row.names(annotation_col)<-colnames(toPlotHmap)

library(RColorBrewer)
pdf("HeatMfuzz_10.pdf",width=8,height = 4)
p1<-pheatmap::pheatmap(toPlotHmap[,order(factorNames,decreasing = T)],scale='row',cluster_rows=F,cluster_cols=F
                       ,color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(15),border_color=NA
                       ,annotation_col = annotation_col, annotation_colors = ann_colors[1])
dev.off()


# pathway enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
Hs.egUniAll <- as.list(org.Hs.egUNIPROT[mappedkeys(org.Hs.egUNIPROT)])
Hs.egUniAllEle<-unlist(Hs.egUniAll)
names(Hs.egUniAllEle)<-Hs.egUniAllEle
m=1;
for(i in 1:length(Hs.egUniAll)){#length(Hs.egUniAll)
  countAdd<-length(Hs.egUniAll[[i]])
  Hs.egUniAllEle[m:(m+countAdd-1)]<-rep(names(Hs.egUniAll)[i],countAdd)
  m<-m+countAdd
}

pdf("mfuzz_t100_2_GO.pdf")
df_UpDw<- data.frame(cl2$cluster,names(cl2$cluster))
colnames(df_UpDw)[1]<-"Group"
df_UpDw$Entrez<-Hs.egUniAllEle[df_UpDw[,2]]
df_UpDwNA0<-df_UpDw[!is.na(df_UpDw$Entrez),]
formula_res2 <- compareCluster(Entrez~Group, data=df_UpDwNA0, fun="enrichGO", OrgDb='org.Hs.eg.db')
dotplot(formula_res2)
dev.off()

pdf("mfuzz_t100_4_GO.pdf")
df_UpDw<- data.frame(cl3$cluster,names(cl3$cluster))
colnames(df_UpDw)[1]<-"Group"
df_UpDw$Entrez<-Hs.egUniAllEle[df_UpDw[,2]]
df_UpDwNA0<-df_UpDw[!is.na(df_UpDw$Entrez),]
formula_res3 <- compareCluster(Entrez~Group, data=df_UpDwNA0, fun="enrichGO", OrgDb='org.Hs.eg.db')
dotplot(formula_res3)
dev.off()

pdf("mfuzz_all_4_GO.pdf")
df_UpDw<- data.frame(cl$cluster,names(cl$cluster))
colnames(df_UpDw)[1]<-"Group"
df_UpDw$Entrez<-Hs.egUniAllEle[df_UpDw[,2]]
df_UpDwNA0<-df_UpDw[!is.na(df_UpDw$Entrez),]
formula_res <- compareCluster(Entrez~Group, data=df_UpDwNA0, fun="enrichGO", OrgDb='org.Hs.eg.db')
dotplot(formula_res)
dev.off()

#UMAP after RF####
library(umap)
pchIDs = seq(10,15)
colIDs = hcl.colors(6, "viridis", rev = TRUE)
names(pchIDs)<-names(colIDs)<-c('N','L','H',"C","PC","CC")

NA_thrs=c(0.4,0.5,0.6)
NA_fill_ratios=c(0.7,0.8,0.9)
NA_thr=0.6
NA_fill_ratio=0.7

pdf("UMAP_NA_thrs0.6fill0.7_RF.pdf",width=42,height = 12)
par(mfcol=c(2,7))
# for(NA_thr in NA_thrs){
  # for(NA_fill_ratio in NA_fill_ratios){
proMatList<-list(proMat,proMat_P1,proMat_P2S1,proMat_P3S1,proMat_P1S1,proMat_P1S2,proMat_P1S3)
names(proMatList)<-c("All","P1","P2","P3","P1S1","P1S2","P1S3")
for (i in 1:length(proMatList)){
  proSel=proMatList[[i]]
  proSel_NA<-proSel[(apply(is.na(proSel),1,sum)/ncol(proSel))<=NA_thr,]
  proSel_NAFill<-proSel_NA
  proSel_NAFill[is.na(proSel_NA)]<-min(proSel_NA,na.rm=T)*NA_fill_ratio
  um<-umap(t(proSel_NAFill))
  subTypesSel<-sapply(colnames(proSel_NAFill),toPlotPar)[3,]
  plot(um$layout[,1],um$layout[,2],xlab='umap1',ylab='umap2'
       ,pch=19,cex=1.5
       ,col=colIDs[subTypesSel]
       ,main = paste0(names(proMatList)[i],"NAthr-",NA_thr,"NAfill-",NA_fill_ratio))
  # text(um$layout[,1],um$layout[,2],col=colIDs[subTypesSel],colnames(proSel_NAFill),cex=0.6)
  
  proSel=proSel[ImpPros,]
  proSel_NA<-proSel[(apply(is.na(proSel),1,sum)/ncol(proSel))<=NA_thr,]
  proSel_NAFill<-proSel_NA
  proSel_NAFill[is.na(proSel_NA)]<-min(proSel_NA,na.rm=T)*NA_fill_ratio
  um<-umap(t(proSel_NAFill))
  subTypesSel<-sapply(colnames(proSel_NAFill),toPlotPar)[3,]
  plot(um$layout[,1],um$layout[,2],xlab='umap1',ylab='umap2'
       ,pch=19,cex=1.5
       ,col=colIDs[subTypesSel]
       ,main = paste0(names(proMatList)[i],"NAthr-",NA_thr,"NAfill-",NA_fill_ratio))
  # text(um$layout[,1],um$layout[,2],col=colIDs[subTypesSel],colnames(proSel_NAFill),cex=0.6)
}
  # }
# }
dev.off()


#centroid velocity
pdf("UMAP_NA_thrs0.6fill0.7_RF_4C.pdf",width=6,height = 6)
proSel=proMatList[[1]]
proSel=proSel[ImpPros,]
proSel_NA<-proSel[(apply(is.na(proSel),1,sum)/ncol(proSel))<=NA_thr,]
proSel_NAFill<-proSel_NA
proSel_NAFill[is.na(proSel_NA)]<-min(proSel_NA,na.rm=T)*NA_fill_ratio
um<-umap(t(proSel_NAFill))
subTypesSel<-sapply(colnames(proSel_NAFill),toPlotPar)[3,]
# pdf("Velocity_map.pdf")

subTypesSel[subTypesSel=="PC"]<-"C"
subTypesSel[subTypesSel=="CC"]<-"C"

plot(um$layout[,1],um$layout[,2],xlab='umap1',ylab='umap2'
     ,pch=19,cex=1.5
     ,col=colIDs[subTypesSel]
     ,main = paste0(names(proMatList)[i],"NAthr-",NA_thr,"NAfill-",NA_fill_ratio))

Median_LocY<-Median_LocX<-rep(NA,4)
names(Median_LocY)<-names(Median_LocX)<-c('N','L','H',"C")
for(i in 1:length(Median_LocY)){
  Median_LocX[i]<-median(um$layout[,1][grepl(names(Median_LocX)[i],subTypesSel)])
  Median_LocY[i]<-median(um$layout[,2][grepl(names(Median_LocY)[i],subTypesSel)])
}
for(i in 1:(length(Median_LocY)-1)){
  arrows(x0=Median_LocX[i],x1=Median_LocX[i+1],
           y0=Median_LocY[i],y1=Median_LocY[i+1],col=colIDs[i],lwd=2)
}
dev.off()

#20230423####
setwd("Figure6")
# Implementation steps
## ● Confirm the positions and orders of punches for five images####
### P1S1####
slideContour<-read.csv("P1S1_raw_contours_slide.csv",stringsAsFactors = F,header = F)
punchContour<-read.csv("P1S1_raw_contours_punch.csv",stringsAsFactors = F,header = F)

{pdf("P1S1_CountourInd.pdf",width=6*2850/800,height=6)
plot(100,100,type='n',xlim=c(0,2850),ylim=c(800,0),axes=F,xlab="",ylab="")
contourGroup<-unique(slideContour$V3)
for(i in 1:length(contourGroup)){
  coutourSel<-slideContour[slideContour$V3==contourGroup[i],]
  if(nrow(coutourSel)>266){
  lines(coutourSel[,1],coutourSel[,2])
  }
}
contourGroup<-unique(punchContour$V3)
punchMeanX<-punchMeanY<-c()
for(i in 1:length(contourGroup)){
  coutourSel<-punchContour[punchContour$V3==contourGroup[i],]
  if(nrow(coutourSel)>27){
    lines(coutourSel[,1],coutourSel[,2],col=i)
    punchMeanX<-c(punchMeanX,mean(coutourSel[,1]))
    punchMeanY<-c(punchMeanY,mean(coutourSel[,2]))
  }
}
# points(punchMeanX,punchMeanY,col=2,pch=19)
text(punchMeanX,punchMeanY,labels = seq(1,length(punchMeanX)))
dev.off()
}
df_punch<-cbind(punchMeanX,punchMeanY)
write.csv(df_punch,"df_punch_P1S1.csv")
df_punch_manual<-read.csv("df_punch_P1S1_manual.csv",stringsAsFactors = F,row.names = 1)
id_all<-sort(unique(df_punch_manual$ID))[-1]
meanX_all<-meanY_all<-c()
for(i in 1:length(id_all)){
  meanX_all[i]<-mean(df_punch_manual$punchMeanX[df_punch_manual$ID==id_all[i]])
  meanY_all[i]<-mean(df_punch_manual$punchMeanY[df_punch_manual$ID==id_all[i]])
}

plot(100,100,type='n',xlim=c(0,2850),ylim=c(800,0),axes=F,xlab="",ylab="")
coutourSel<-slideContour[slideContour$V3==max(slideContour$V3),1:2]
lines(coutourSel[,1],coutourSel[,2])
points(meanX_all,meanY_all,cex=2)
write.csv(coutourSel,"P1S1_2850-800_slideContour.csv",row.names = F)
df_meanLoc<-cbind(meanX_all,meanY_all)
row.names(df_meanLoc)<-id_all
write.csv(df_meanLoc,"P1S1_2850-800_punchLoc.csv")
#P1S2
slideContour<-read.csv("P1S2_raw_contours_slide.csv",stringsAsFactors = F,header = F)
punchContour<-read.csv("P1S2_raw_contours_punch.csv",stringsAsFactors = F,header = F)

{pdf("P1S2_CountourInd.pdf",width=6*2866/871,height=6)
  plot(100,100,type='n',xlim=c(0,2866),ylim=c(871,0),axes=F,xlab="",ylab="")
  contourGroup<-unique(slideContour$V3)
  coutourSel<-slideContour[slideContour$V3==(which.max(table(slideContour$V3))-1),]
  lines(coutourSel[,1],coutourSel[,2])
  
  contourGroup<-unique(punchContour$V3)
  punchMeanX<-punchMeanY<-c()
  for(i in 1:length(contourGroup)){
    coutourSel<-punchContour[punchContour$V3==contourGroup[i],]
    if(nrow(coutourSel)>25){
      lines(coutourSel[,1],coutourSel[,2],col=i)
      punchMeanX<-c(punchMeanX,mean(coutourSel[,1]))
      punchMeanY<-c(punchMeanY,mean(coutourSel[,2]))
    }
  }
  text(punchMeanX,punchMeanY,labels = seq(1,length(punchMeanX)))
  dev.off()
}
df_punch<-cbind(punchMeanX,punchMeanY)
write.csv(df_punch,"df_punch_P1S2.csv")

df_punch_manual<-read.csv("df_punch_P1S2_manual.csv",stringsAsFactors = F,row.names = 1)
id_all<-sort(unique(df_punch_manual$ID))[-1]
meanX_all<-meanY_all<-c()
for(i in 1:length(id_all)){
  meanX_all[i]<-mean(df_punch_manual$punchMeanX[df_punch_manual$ID==id_all[i]])
  meanY_all[i]<-mean(df_punch_manual$punchMeanY[df_punch_manual$ID==id_all[i]])
}

plot(100,100,type='n',xlim=c(0,2866),ylim=c(871,0),axes=F,xlab="",ylab="")
coutourSel<-slideContour[slideContour$V3==(which.max(table(slideContour$V3))-1),]
lines(coutourSel[,1],coutourSel[,2])
points(meanX_all,meanY_all,cex=2)
write.csv(coutourSel,"P1S2_2866-871_slideContour.csv",row.names = F)
df_meanLoc<-cbind(meanX_all,meanY_all)
row.names(df_meanLoc)<-id_all
write.csv(df_meanLoc,"P1S2_2866-871_punchLoc.csv")

#P1S3
slideContour<-read.csv("P1S3_raw_contours_slide.csv",stringsAsFactors = F,header = F)
punchContour<-read.csv("P1S3_raw_contours_punch.csv",stringsAsFactors = F,header = F)
imgWidth=2697;imgHeight=841

{pdf("P1S3_CountourInd.pdf",width=6*imgWidth/imgHeight,height=6)
  plot(100,100,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab="")
  contourGroup<-unique(slideContour$V3)
  coutourSel<-slideContour[slideContour$V3==(which.max(table(slideContour$V3))-1),]
  lines(coutourSel[,1],coutourSel[,2])
  
  contourGroup<-unique(punchContour$V3)
  punchMeanX<-punchMeanY<-c()
  for(i in 1:length(contourGroup)){
    coutourSel<-punchContour[punchContour$V3==contourGroup[i],]
    if(nrow(coutourSel)>25){
      lines(coutourSel[,1],coutourSel[,2],col=i)
      punchMeanX<-c(punchMeanX,mean(coutourSel[,1]))
      punchMeanY<-c(punchMeanY,mean(coutourSel[,2]))
    }
  }
  text(punchMeanX,punchMeanY,labels = seq(1,length(punchMeanX)))
  dev.off()
}
df_punch<-cbind(punchMeanX,punchMeanY)
write.csv(df_punch,"df_punch_P1S3.csv")

df_punch_manual<-read.csv("df_punch_P1S3_manual.csv",stringsAsFactors = F,row.names = 1)
id_all<-sort(unique(df_punch_manual$ID))[-1]
meanX_all<-meanY_all<-c()
for(i in 1:length(id_all)){
  meanX_all[i]<-mean(df_punch_manual$punchMeanX[df_punch_manual$ID==id_all[i]])
  meanY_all[i]<-mean(df_punch_manual$punchMeanY[df_punch_manual$ID==id_all[i]])
}

plot(100,100,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab="")
coutourSel<-slideContour[slideContour$V3==(which.max(table(slideContour$V3))-1),]
lines(coutourSel[,1],coutourSel[,2])
points(meanX_all,meanY_all,cex=1.75)
write.csv(coutourSel,paste0("P1S3_",imgWidth,"-",imgHeight,"_slideContour.csv"),row.names = F)
df_meanLoc<-cbind(meanX_all,meanY_all)
row.names(df_meanLoc)<-id_all
write.csv(df_meanLoc,paste0("P1S3_",imgWidth,"-",imgHeight,"_punchLoc.csv"))

#P2S1
slideContour<-read.csv("P2S1_raw_contours_slide.csv",stringsAsFactors = F,header = F)
punchContour<-read.csv("P2S1_raw_contours_punch.csv",stringsAsFactors = F,header = F)
imgWidth=2939;imgHeight=3609

{pdf("P2S1_CountourInd.pdf",width=6*imgHeight/imgWidth,height=6)
  par(mar=c(0,0,0,0))
  plot(0,0,type='n',xlim=c(0,imgHeight),ylim=c(imgWidth,0),axes=F,xlab="",ylab="")
  coutourSel<-slideContour[slideContour$V3==(which.max(table(slideContour$V3))-1),]
  lines(coutourSel[,1],coutourSel[,2])
  
  contourGroup<-unique(punchContour$V3)
  punchMeanX<-punchMeanY<-c()
  for(i in 1:length(contourGroup)){
    coutourSel<-punchContour[punchContour$V3==contourGroup[i],]
    if(nrow(coutourSel)>25){
      lines(coutourSel[,1],coutourSel[,2],col=i)
      punchMeanX<-c(punchMeanX,mean(coutourSel[,1]))
      punchMeanY<-c(punchMeanY,mean(coutourSel[,2]))
    }
  }
  text(punchMeanX,punchMeanY,labels = seq(1,length(punchMeanX)),cex=0.3)
  dev.off()
}
df_punch<-cbind(punchMeanX,punchMeanY)
write.csv(df_punch,"df_punch_P2S1.csv")

df_punch_manual<-read.csv("df_punch_P2S1_manual.csv",stringsAsFactors = F,row.names = 1)
id_all<-sort(unique(df_punch_manual$ID))[-1]
meanX_all<-meanY_all<-c()
for(i in 1:length(id_all)){
  meanX_all[i]<-mean(df_punch_manual$punchMeanX[df_punch_manual$ID==id_all[i]])
  meanY_all[i]<-mean(df_punch_manual$punchMeanY[df_punch_manual$ID==id_all[i]])
}

plot(100,100,type='n',xlim=c(0,imgHeight),ylim=c(imgWidth,0),axes=F,xlab="",ylab="")
coutourSel<-slideContour[slideContour$V3==(which.max(table(slideContour$V3))-1),]
lines(coutourSel[,1],coutourSel[,2])
points(meanX_all,meanY_all,cex=1.75)
write.csv(coutourSel,paste0("P2S1_",imgHeight,"-",imgWidth,"_slideContour.csv"),row.names = F)
df_meanLoc<-cbind(meanX_all,meanY_all)
row.names(df_meanLoc)<-id_all
write.csv(df_meanLoc,paste0("P2S1_",imgHeight,"-",imgWidth,"_punchLoc.csv"))

#P3S1
slideContour<-read.csv("P3S1_raw_contours_slide.csv",stringsAsFactors = F,header = F)
punchContour<-read.csv("P3S1_raw_contours_punch.csv",stringsAsFactors = F,header = F)
imgWidth=3954;imgHeight=3405

{pdf("P3S1_CountourInd.pdf",width=6*imgWidth/imgHeight,height=6)
  par(mar=c(0,0,0,0))
  plot(0,0,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab="")
  coutourSel<-slideContour[slideContour$V3==(which.max(table(slideContour$V3))-1),]
  lines(coutourSel[,1],coutourSel[,2])
  
  contourGroup<-unique(punchContour$V3)
  punchMeanX<-punchMeanY<-c()
  for(i in 1:length(contourGroup)){
    coutourSel<-punchContour[punchContour$V3==contourGroup[i],]
    if(nrow(coutourSel)>15){
      lines(coutourSel[,1],coutourSel[,2],col=i)
      punchMeanX<-c(punchMeanX,mean(coutourSel[,1]))
      punchMeanY<-c(punchMeanY,mean(coutourSel[,2]))
    }
  }
  text(punchMeanX,punchMeanY,labels = seq(1,length(punchMeanX)))
  dev.off()
}
df_punch<-cbind(punchMeanX,punchMeanY)
write.csv(df_punch,"df_punch_P3S1.csv")

df_punch_manual<-read.csv("df_punch_P3S1_manual.csv",stringsAsFactors = F,row.names = 1)
id_all<-sort(unique(df_punch_manual$ID))[-1]
meanX_all<-meanY_all<-c()
for(i in 1:length(id_all)){
  meanX_all[i]<-mean(df_punch_manual$punchMeanX[df_punch_manual$ID==id_all[i]])
  meanY_all[i]<-mean(df_punch_manual$punchMeanY[df_punch_manual$ID==id_all[i]])
}

plot(100,100,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab="")
coutourSel<-slideContour[slideContour$V3==(which.max(table(slideContour$V3))-1),]
lines(coutourSel[,1],coutourSel[,2])
points(meanX_all,meanY_all,cex=1.75)
write.csv(coutourSel,paste0("P3S1_",imgWidth,"-",imgHeight,"_slideContour.csv"),row.names = F)
df_meanLoc<-cbind(meanX_all,meanY_all)
row.names(df_meanLoc)<-id_all
write.csv(df_meanLoc,paste0("P3S1_",imgWidth,"-",imgHeight,"_punchLoc.csv"))

# ● Using whole proteins for an unsupervised clustering to separate regions for each slide.
slideContourFiles<-list.files(pattern="_slideContour.csv")
punchLocFiles<-list.files(pattern="_punchLoc.csv")
names(slideContourFiles)<-names(punchLocFiles)<-c("P1S1","P1S2","P1S3","P2S1","P3S1")
pchIDs = seq(10,15)
colIDs = hcl.colors(6, "viridis", rev = TRUE)
names(pchIDs)<-names(colIDs)<-c('N','L','H',"C","PC","CC")
#P1S1
proMat_P1S1_NA[is.na(proMat_P1S1)]
# for(NA_thr in NA_thrs){
  # for(NA_fill_ratio in NA_fill_ratios){
NA_thr=NA_fill_ratio=0.8
proMatList<-list(proMat_P1S1,proMat_P1S2,proMat_P1S3,proMat_P2S1,proMat_P3S1)#proMat,proMat_P1,
names(proMatList)<-c("P1S1","P1S2","P1S3","P2S1","P3S1")
for (i in 1:3){#length(proMatList)
  proSel=proMatList[[i]]
  proSel_NA<-proSel[(apply(is.na(proSel),1,sum)/ncol(proSel))<=NA_thr,]
  proSel_NAFill<-proSel_NA
  proSel_NAFill[is.na(proSel_NA)]<-min(proSel_NA,na.rm=T)*NA_fill_ratio
  proSel_NAFillZ<-proSel_NAFill
  for(m in 1:nrow(proSel_NAFill)){
    proSel_NAFillZ[m,]<-scale(as.numeric(proSel_NAFill[m,]))[,1]
  }
  
  library(cluster)
  set.seed(123)
  cenNum<-length(unique(sapply(colnames(proSel_NAFillZ),toPlotPar)[3,]))
  kmeans_model <- kmeans(t(proSel_NAFillZ), centers = cenNum)
  dist_Punch<-as.matrix(dist(t(proSel_NAFillZ)))
  dist_Punch_Sel<-dist_Punch
  cutNet<-100
  dist_Punch_Sel[dist_Punch>cutNet]=0
  dist_Punch_Sel[dist_Punch<=cutNet]=1
  
  punchLoc<-read.csv(punchLocFiles[i],stringsAsFactors = F,row.names = 1)
  slideContour<-read.csv(slideContourFiles[i],stringsAsFactors = F)
  imgHeight=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[2])
  imgWidth=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[1])
  # pdf(paste0("Cluster_Ellipse0.95_",names(proMatList)[i],".pdf"),width=6*imgWidth/imgHeight,height=6)
  par(mar=c(1,1,4,1))
  plot(0,0,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab=""
       ,main=names(proMatList)[i])
  lines(slideContour[,1],slideContour[,2])
  for(n in 1:(nrow(dist_Punch_Sel)-1)){
    for(m in (n+1):ncol(dist_Punch_Sel)){
      if(dist_Punch_Sel[m,n]==1){
        lines(punchLoc[c(row.names(dist_Punch_Sel)[m],colnames(dist_Punch_Sel)[n]),1],
              punchLoc[c(row.names(dist_Punch_Sel)[m],colnames(dist_Punch_Sel)[n]),2])
      }
    }
  }
  
  points(punchLoc[,1],punchLoc[,2],pch=kmeans_model$cluster[row.names(punchLoc)]+20
         ,cex=2,bg=colIDs[sapply(row.names(punchLoc),toPlotPar)[3,]])
  # centroidX<-centroidY<-rep(0,cenNum)
  for(n in 1:cenNum){
    punchSel<-punchLoc[names(kmeans_model$cluster)[kmeans_model$cluster==n],]
    # polygon(punchSel[,1],punchSel[,2])
    library(car)
    if(nrow(punchSel)>1){
      d1=dataEllipse(punchSel$meanX_all,punchSel$meanY_all,draw=F,levels=0.75)
      lines(d1)
    }
  }

  # dev.off()
}

#  Using dys-regulated protein visualizations as networks for different slides####

proMatList<-list(proMat_P1S1,proMat_P1S2,proMat_P1S3,proMat_P2S1,proMat_P3S1)#proMat,proMat_P1,
names(proMatList)<-c("P1S1","P1S2","P1S3","P2S1","P3S1")
for (i in 1:length(proMatList)){
  proSel=proMatList[[i]]
  proSelZ<-proSel
  for(m in 1:nrow(proSel)){
    proSelZ[m,]<-scale(as.numeric(proSel[m,]))[,1]
  }
  punchLoc<-read.csv(punchLocFiles[i],stringsAsFactors = F,row.names = 1)
  slideContour<-read.csv(slideContourFiles[i],stringsAsFactors = F)
  imgHeight=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[2])
  imgWidth=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[1])
  
  # pdf(paste0("HeatAllZ_",names(proMatList)[i],".pdf"),width=120,height=100)
  # par(mfrow=c(40,30))
  for(n in 1:nrow(proSelZ)){#nrow(proSelZ)
    par(mar=c(0,0,3,0))
    plot(0,0,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab=""
         ,main=row.names(proSelZ)[n])
    lines(slideContour[,1],slideContour[,2])
    pchAll<-rep(18,ncol(proSelZ))
    pchAll[is.na(proSelZ[n,])]<-23
    colAll<-colGenMax(proSelZ[n,],hclPatl)
    colAll[is.na(colAll)]<-1
    points(punchLoc[,1],punchLoc[,2],cex=2.5,pch=pchAll,col=colAll)
    print(n)
  }
  rgl.postscript('3dplot.pdf', fmt = 'pdf')  
}

# ● For the Patient 1 can draw a 3D image.
library(rgl)
proMatList<-list(proMat_P1S1,proMat_P1S2,proMat_P1S3)
names(proMatList)<-c("P1S1","P1S2","P1S3")
punchLocs<-contourLines<-list()
for(i in 1:length(proMatList)){
  punchLocs[[i]]<-read.csv(punchLocFiles[i],stringsAsFactors = F,row.names = 1)
  contourLines[[i]]<-read.csv(slideContourFiles[i],stringsAsFactors = F)
}
# for(m in 1:nrow(proMat)){}
m = which(row.names(proMat)=="P01116")
plot3d(0,0,0,type='n', decorate=F,zlim=c(0,1500),xlim=c(0,2900),ylim=c(0,870)
       ,main=row.names(proMat)[m],aspect=10)
for (i in 1:length(proMatList)){
  proSel=proMatList[[i]]
  proSelZ<-scale(as.numeric(proSel[m,]))[,1]
  colAll<-colGenMax(proSelZ,heatmap_colors)
  colAll[is.na(colAll)]<-NA
  
  plot3d(punchLocs[[i]][,1], punchLocs[[i]][,2], i*500, type="s", col=colAll, size=1.5,add=T)
  plot3d(contourLines[[i]][,1],contourLines[[i]][,2], i*500, type="l",add=T)
  
}
view3d(theta=-5, phi=-60, fov=0)
rgl.postscript(paste0("Heat3D_",row.names(proMat)[m],".pdf"), fmt = 'pdf')


#HE partition####
# for(i in 1:length(punchLocFiles)){
#   punchLoc<-read.csv(punchLocFiles[i],stringsAsFactors = F,row.names = 1)
#   slideContour<-read.csv(slideContourFiles[i],stringsAsFactors = F)
#   imgHeight=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[2])
#   imgWidth=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[1])
#   pdf(paste0(names(punchLocFiles)[i],"_OutCountour.pdf"),width=6*imgWidth/imgHeight,height=6)
#   par(mar=c(0,0,0,0))
#   plot(0,0,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab="")
#   lines(slideContour[,1],slideContour[,2])
#   # pchAll<-rep(18,ncol(proSelZ))
#   # pchAll[is.na(proSelZ[n,])]<-23
#   # colAll<-colGenMax(proSelZ[n,],hclPatl)
#   # colAll[is.na(colAll)]<-1
#   # points(punchLoc[,1],punchLoc[,2],cex=2.5)
#   dev.off()
# }
contoursHE<-read.csv("P1S1_contours_region.csv",stringsAsFactors = F,header = F)
contoursHE_Sel<-contoursHE[contoursHE$V3%%2==0,]
contoursHE_Sel_P1S3<-contoursHE_Sel_P1S2<-contoursHE_Sel_P1S1<-contoursHE_Sel
par(mar=c(0,0,0,0))
# for(i in 3:length(punchLocFiles[1:3])){
#   punchLoc<-read.csv(punchLocFiles[i],stringsAsFactors = F,row.names = 1)
#   slideContour<-read.csv(slideContourFiles[i],stringsAsFactors = F)
#   imgHeight=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[2])
#   imgWidth=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[1])
#   par(mar=c(0,0,0,0))
#   plot(0,0,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab="")
#   polygon(slideContour[,1],slideContour[,2],border ="grey")
#   points(punchLoc[,1],punchLoc[,2],cex=2.5)
  # col_P1S1<-colIDs[c("L","H","L","N","PC","CC")]
  # col_P1S2<-colIDs[c("L","H","L","N","PC","CC")]
  # col_P1S3<-colIDs[c("L","H","L","N","PC","CC")]
  # for(n in c(1:3,5:6)){
    # lines(contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,1]/2.02-246,
    #       contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,2]/2.25+70)
    # contoursHE_Sel_P1S3[contoursHE_Sel$V3==(n-1)*2,1]<-contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,1]/2.02-246
    # contoursHE_Sel_P1S3[contoursHE_Sel$V3==(n-1)*2,2]<-contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,2]/2.25+70
  # }
  # contoursHE_Sel_P1S3[contoursHE_Sel$V3==(4-1)*2,1]<-contoursHE_Sel[contoursHE_Sel$V3==(4-1)*2,1]/2.02-296
  # contoursHE_Sel_P1S3[contoursHE_Sel$V3==(4-1)*2,2]<-contoursHE_Sel[contoursHE_Sel$V3==(4-1)*2,2]/2.1+20
  # lines(contoursHE_Sel[contoursHE_Sel$V3==(4-1)*2,1]/2.02-296,
  #       contoursHE_Sel[contoursHE_Sel$V3==(4-1)*2,2]/2.1+20)
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==0,1]<-contoursHE_Sel[contoursHE_Sel$V3==0,1]/2.1+10
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==2,1]<-contoursHE_Sel[contoursHE_Sel$V3==2,1]/2.1-40
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==4,1]<-contoursHE_Sel[contoursHE_Sel$V3==4,1]/2.1-40
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==6,1]<-contoursHE_Sel[contoursHE_Sel$V3==6,1]/2.1-50
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==8,1]<-contoursHE_Sel[contoursHE_Sel$V3==8,1]/2.1-50
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==10,1]<-contoursHE_Sel[contoursHE_Sel$V3==10,1]/2-176
  # 
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==0,2]<-contoursHE_Sel[contoursHE_Sel$V3==0,2]/2.11+90
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==2,2]<-contoursHE_Sel[contoursHE_Sel$V3==2,2]/2.11+20
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==4,2]<-contoursHE_Sel[contoursHE_Sel$V3==4,2]/2.11+20
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==6,2]<-contoursHE_Sel[contoursHE_Sel$V3==6,2]/2.01+70
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==8,2]<-contoursHE_Sel[contoursHE_Sel$V3==8,2]/2.11+40
  # contoursHE_Sel_P1S2[contoursHE_Sel$V3==10,2]<-contoursHE_Sel[contoursHE_Sel$V3==10,2]/2.27+60
  
  # lines(contoursHE_Sel[contoursHE_Sel$V3==(1-1)*2,1]/2.1+10,contoursHE_Sel[contoursHE_Sel$V3==(1-1)*2,2]/2.11+90,col=1)
  # lines(contoursHE_Sel[contoursHE_Sel$V3==(2-1)*2,1]/2.1-40,contoursHE_Sel[contoursHE_Sel$V3==(2-1)*2,2]/2.11+20,col=2)
  # lines(contoursHE_Sel[contoursHE_Sel$V3==(3-1)*2,1]/2.1-40,contoursHE_Sel[contoursHE_Sel$V3==(3-1)*2,2]/2.11+20,col=3)
  # lines(contoursHE_Sel[contoursHE_Sel$V3==(4-1)*2,1]/2.1-50,contoursHE_Sel[contoursHE_Sel$V3==(4-1)*2,2]/2.01+70,col=4)
  # lines(contoursHE_Sel[contoursHE_Sel$V3==(5-1)*2,1]/2.1-50,contoursHE_Sel[contoursHE_Sel$V3==(5-1)*2,2]/2.11+40,col=5)
  # lines(contoursHE_Sel[contoursHE_Sel$V3==(6-1)*2,1]/2-176,contoursHE_Sel[contoursHE_Sel$V3==(6-1)*2,2]/2.27+60,col=6)

  # for(n in 1:length(unique(contoursHE_Sel$V3))){
    # lines(contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,1]/2.02-246,
    #       contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,2]/2.25+70)
    # contoursHE_Sel_P1S1[contoursHE_Sel$V3==(n-1)*2,1]<-contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,1]/2.02-206
    # contoursHE_Sel_P1S1[contoursHE_Sel$V3==(n-1)*2,2]<-contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,2]/2.25
    # polygon(contoursHE_Sel_P1S1[contoursHE_Sel_P1S1$V3==(n-1)*2,1],
          # contoursHE_Sel_P1S1[contoursHE_Sel_P1S2$V3==(n-1)*2,2],border=col_P1S1[n])
    # polygon(contoursHE_Sel_P1S2[contoursHE_Sel_P1S2$V3==(n-1)*2,1],
    #       contoursHE_Sel_P1S2[contoursHE_Sel_P1S2$V3==(n-1)*2,2],border=col_P1S2[n])
    # polygon(contoursHE_Sel_P1S3[contoursHE_Sel_P1S3$V3==(n-1)*2,1],
    #       contoursHE_Sel_P1S3[contoursHE_Sel_P1S3$V3==(n-1)*2,2],border=col_P1S3[n])
  # }
# }
write.csv(contoursHE_Sel_P1S1,"contoursHE_Sel_P1S1.csv",row.names = F)
write.csv(contoursHE_Sel_P1S2,"contoursHE_Sel_P1S2.csv",row.names = F)
write.csv(contoursHE_Sel_P1S3,"contoursHE_Sel_P1S3.csv",row.names = F)
contoursHEs<-list(contoursHE_Sel_P1S1,contoursHE_Sel_P1S2,contoursHE_Sel_P1S3)
names(contoursHEs)<-c("P1S1","P1S2","P1S3")

# re-draw the 2D visualizations with outlines
m = which(row.names(proMat)=="P01116")
# re-draw the 2D visualizations
# for(n in 1:nrow(proSelZ)){#nrow(proSelZ)
colP1<-c("#7ED35750","#00B28B40","#7ED35750","#FDE333A0","#27498380","#4B005580" )
# colP1<-c("#C7E9C080","#74C47680","#C7E9C080","#EDF8E980","#31A35480","#006D2C80" )

# brewer.pal(6,"Greens")
# "#EDF8E9" "#C7E9C0" "#A1D99B" "#74C476" "#31A354" "#006D2C"
par(mar=c(0,0,0,0),mfrow=c(3,1))
for (i in 1:length(proMatList[1:3])){
  proSel=proMatList[[i]]
  proSelZ<-scale(as.numeric(proSel[m,]))[,1]
  colAll<-colGenMax(proSelZ,hclPatl)
  colAll[is.na(colAll)]<-NA
  imgHeight=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[2])
  imgWidth=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[1])
  plot(0,0,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab=""
       ,main="")#row.names(proMat)[m]
  polygon(contourLines[[i]][,1],contourLines[[i]][,2],col="#D1ACC580",border=NA)
  pchAll<-rep(21,length(proSelZ))
  pchAll[is.na(proSelZ)]<-23
  colAll<-colGenMax(proSelZ,heatmap_colors)
  colAll[is.na(colAll)]<-1
  contoursHE_Sel<-contoursHEs[[i]]
  for(n in 1:length(unique(contoursHE_Sel$V3))){
    polygon(contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,1],
            contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,2],col=colP1[n],border=colP1[n])
  }
  points(punchLocs[[i]][,1],punchLocs[[i]][,2],cex=3.5,pch=pchAll,bg=colAll)
  
}
#three stacks####
imgHeights<-imgWidths<-c()
for(i in 1:length(slideContourFiles)){
  imgHeights[i]=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[2])
  imgWidths[i]=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[i],"_"))[2],"-"))[1])
}
proMats<-list(proMat_P1S1,proMat_P1S2,proMat_P1S3,proMat_P2S1,proMat_P3S1)
slideContourFiles<-list.files("../Figure6",pattern="_slideContour.csv")
punchLocFiles<-list.files("../Figure6",pattern="_punchLoc.csv")
punchLocs<-slideContours<-list()
for(i in 1:5){
  punchLocs[[i]]<-read.csv(paste0("../Figure6/",punchLocFiles[i]),stringsAsFactors = F,row.names = 1)
  slideContours[[i]]<-read.csv(paste0("../Figure6/",slideContourFiles[i]),stringsAsFactors = F)
}
load("../ProtomExAct/contoursHEs.Rdata")

SelTop4<-c("P31946", "P43490","O43865","P85298")
for(m in 1:length(SelTop4)){
  pdf(paste0("Three_stacks_",SelTop4[m],"_",GnMatch[SelTop4[m],"Gene.Names.Main"],".pdf"),width=9,height=6)
  par(mar=c(0,0,0,0),mfcol=c(3,1))
  for(i in 1:3){#P1
    proSel<-SelTop4[m]#"P01116"
    proSelZ<-scale(as.numeric(proMats[[i]][proSel,]))[,1]
    names(proSelZ)<-colnames(proMats[[i]])
    
    punchLoc<-punchLocs[[i]]
    
    colP1<-c("#7ED35750","#00B28B40","#7ED35750","#FDE333A0","#27498380","#4B005580" )
    
    plot(0,0,type='n',xlim=c(0,imgWidths[i]),ylim=c(imgHeights[i],0),axes=F,xlab="",ylab="",main=paste(SelTop4[n], GnMatch[SelTop4[m],"Gene.Names.Main"],sep="_"))
    polygon(slideContours[[i]][,1],slideContours[[i]][,2],col=NA,border=1)#"#9195BA80",
    pchAll<-rep(21,length(proSelZ))
    pchAll[is.na(proSelZ)]<-21
    colAll<-colGenMax(proSelZ,rev(brewer.pal(7, "RdYlBu")))
    colAll[is.na(colAll)]<-NA
    contoursHE_Sel<-contoursHEs[[i]]
    for(n in 1:length(unique(contoursHE_Sel$V3))){
      polygon(contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,1],
              contoursHE_Sel[contoursHE_Sel$V3==(n-1)*2,2],border=colP1[n],col=NA, lwd = 4)#col=colP1[n],
    }#GnMatch
    points(punchLocs[[i]][,1],punchLocs[[i]][,2],cex=3.5,pch=19,col="#9195BA80",border =NA)
    points(punchLocs[[i]][,1],punchLocs[[i]][,2],cex=2,pch=21,bg=colAll,border =NA)
    
  }
  dev.off()
}


# color_range <- colorRampPalette(c("blue", "white", "red"))
# heatmap_colors <-colorRampPalette(c("blue", "white", "red"))(9)
pdf("ThreeStackColorBar.pdf")
heatmap_colors<-rev(brewer.pal(7, "RdYlBu"))
colBar<-rep(1,length(heatmap_colors))
names(colBar)<-c(-1.8,-1.2,-0.6,0,0.6,1.2,1.8)
barplot(colBar, col = heatmap_colors,space=0,border=NA)
dev.off()

pdf("KEGGColorBar.pdf")
heatmap_colors1<-rev(brewer.pal(9, "RdYlBu"))
colBar<-rep(1,length(heatmap_colors1))
names(colBar)<-seq(-1,1,length.out=9)
barplot(colBar, col = heatmap_colors1,space=0,border=NA)
dev.off()

colGenMax<-function(vectorIn,ColIn){
  maxValue<-1.8
  vectorIn[which(vectorIn>maxValue)]<-maxValue
  vectorIn[which(vectorIn< -maxValue)]<- -maxValue
  
  diffValue=maxValue*2/(length(ColIn)-1)
  intIndex = as.integer((vectorIn+maxValue)/diffValue)+1
  return(ColIn[intIndex])
}

colGenMax1<-function(vectorIn,ColIn){
  maxValue<-1
  vectorIn[which(vectorIn>maxValue)]<-maxValue
  vectorIn[which(vectorIn< -maxValue)]<- -maxValue
  
  diffValue=maxValue*2/(length(ColIn)-1)
  intIndex = as.integer((vectorIn+maxValue)/diffValue)+1
  return(ColIn[intIndex])
}



#3D####
# m = which(row.names(proMat)=="P01116")
{
  plot3d(0,0,0,type='n', decorate=F,zlim=c(0,1500),xlim=c(0,2900),ylim=c(0,870*2)
       ,main=row.names(proMat)[m],aspect=10)
# rgl.postscript(paste0("Heat3D_",row.names(proMat)[m],".pdf"), fmt = 'pdf')
library(cluster)
set.seed(123)
# proSel=proMat_P1
# proSel_NA<-proSel[(apply(is.na(proSel),1,sum)/ncol(proSel))<=NA_thr,]
# proSel_NAFill<-proSel_NA
# proSel_NAFill[is.na(proSel_NA)]<-min(proSel_NA,na.rm=T)*NA_fill_ratio
# proSel_NAFillZ<-proSel_NAFill
# for(m in 1:nrow(proSel_NAFill)){
#   proSel_NAFillZ[m,]<-scale(as.numeric(proSel_NAFill[m,]))[,1]
# }
cenNum<-4#length(unique(sapply(colnames(proSel_NAFillZ),toPlotPar)[3,]))
kmeans_model <- kmeans(t(proSel_NAFillZ), centers = cenNum)
# dist_Punch<-as.matrix(dist(t(proSel_NAFillZ)))
# dist_Punch_Sel<-dist_Punch
# cutNet<-80
# dist_Punch_Sel[dist_Punch>cutNet]=0
# dist_Punch_Sel[dist_Punch<=cutNet]=1
# for(n in 1:(nrow(dist_Punch_Sel)-1)){
#   for(m in (n+1):ncol(dist_Punch_Sel)){
#     if(dist_Punch_Sel[m,n]==1){
#       x1 = punchLocs[[as.numeric(unlist(strsplit(row.names(dist_Punch_Sel)[m],""))[4])]][row.names(dist_Punch_Sel)[m],1]
#       x2 = punchLocs[[as.numeric(unlist(strsplit(colnames(dist_Punch_Sel)[n],""))[4])]][colnames(dist_Punch_Sel)[n],1]
#       y1 = punchLocs[[as.numeric(unlist(strsplit(row.names(dist_Punch_Sel)[m],""))[4])]][row.names(dist_Punch_Sel)[m],2]
#       y2 = punchLocs[[as.numeric(unlist(strsplit(colnames(dist_Punch_Sel)[n],""))[4])]][colnames(dist_Punch_Sel)[n],2]
#       z1 = as.numeric(unlist(strsplit(row.names(dist_Punch_Sel)[m],""))[4])*500
#       z2 = as.numeric(unlist(strsplit(colnames(dist_Punch_Sel)[n],""))[4])*500
#       plot3d(c(x1,x2),c(y1,y2),c(z1,z2), type="l",add=T,col="#D1ACC580")
#     }
#   }
# }
punchLocAll<-rbind(punchLocs[[1]],punchLocs[[2]],punchLocs[[3]])
punchLocAll$Z<-c(rep(1,nrow(punchLocs[[1]])),rep(2,nrow(punchLocs[[2]])),rep(3,nrow(punchLocs[[2]])))
punchLocAll$Cluster<-kmeans_model$cluster[row.names(punchLocAll)]
clusterX<-clusterY<-clusterZ<-c()
for(i in 1:cenNum){
  clusterX[i]<-median(punchLocAll[punchLocAll$Cluster==i,1])
  clusterY[i]<-median(punchLocAll[punchLocAll$Cluster==i,2])*2
  clusterZ[i]<-median(punchLocAll[punchLocAll$Cluster==i,3])*500
}
for (i in 1:length(punchLocs)){
  # proSel=proMatList[[i]]
  # proSelZ<-scale(as.numeric(proSel[m,]))[,1]
  # colAll<-colGenMax(proSelZ,heatmap_colors)
  # colAll[is.na(colAll)]<-NA
  # kmeans_model$cluster[row.names(punchLocs[[i]])]
  plot3d(contourLines[[i]][,1],contourLines[[i]][,2]*2, i*500, type="l",add=T,col="orange")
  # polygon3d(contourLines[[i]][,1],contourLines[[i]][,2]*2, i*500, type="l",add=T,col="#9195BA80")
  # polygon3d(contourLines[[i]][,1],contourLines[[i]][,2]*2,rep(i*500,nrow(contourLines[[i]])),col="#9195BA80")
  # plot3d(punchLocs[[i]][,1], punchLocs[[i]][,2], i*500, type="s", col=kmeans_model$cluster[row.names(punchLocs[[i]])], size=1.5,add=T)
  pch3d(punchLocs[[i]][,1], punchLocs[[i]][,2]*2, rep(i*500,nrow(punchLocs[[i]])),pch=kmeans_model$cluster[row.names(punchLocs[[i]])]+20, cex=0.25,add=T,bg = colIDs[sapply(row.names(punchLocs[[i]]),toPlotPar)[3,]])
}
arrow3d(p0=c(clusterX[3],clusterY[3],clusterZ[3]),p1=c(clusterX[2],clusterY[2],clusterZ[2]), type = "flat",col="black",barblen =0.05,thickness=123)
arrow3d(p0=c(clusterX[2],clusterY[2],clusterZ[2]),p1=c(clusterX[4],clusterY[4],clusterZ[4]), type = "flat",col="black",barblen =0.05,thickness=123)
arrow3d(p0=c(clusterX[4],clusterY[4],clusterZ[4]),p1=c(clusterX[1],clusterY[1],clusterZ[1]), type = "flat",col="black",barblen =0.05,thickness=123)
# arrow3d(p0=c(clusterX[4],clusterY[4],clusterZ[4]),p1=c(clusterX[5],clusterY[5],clusterZ[5]), type = "lines",col="black",barblen =0.01)
view3d(theta=-5, phi=-60, fov=0)
}
rgl.postscript("3D_20230701v2.pdf",fmt="pdf")


# # centroidX<-centroidY<-rep(0,cenNum)
# for(n in 1:cenNum){
#   punchSel<-punchLoc[names(kmeans_model$cluster)[kmeans_model$cluster==n],]
#   # polygon(punchSel[,1],punchSel[,2])
#   library(car)
#   if(nrow(punchSel)>1){
#     d1=dataEllipse(punchSel$meanX_all,punchSel$meanY_all,draw=F,levels=0.75)
#     lines(d1)
#   }
# }

#  Visualize with z scores for each image as 2D visualization####.
# count_3 <- c(20, 50, 30)
# pie(count_3, labels = paste0(count_3, "%"))
# pie(count_3, labels = paste0(count_3, "%"),add=T,radius=0.5)
# plot(0, xlim = c(-1, 1), ylim = c(-1, 1), type = "n", xlab = "", ylab = "", asp = 1)
# symbols(0, 0, circles = 0.5, add = TRUE, inches = FALSE, bg = "blue")
# symbols(0, 0, circles = 0.3, add = TRUE, inches = FALSE, bg = "red")

# data <- c(25, 20, 15, 10, 30)
# labels <- c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5")
# angles <- 360 * data / sum(data)
# radius <- sqrt(data / sum(data))
# plot.new()
# plot.window(c(-1, 1), c(-1, 1))
# # axis(1, at = seq(-1, 1, by = 0.2), labels = NA, lwd = 0, lwd.ticks = 1)
# # axis(2, at = seq(-1, 1, by = 0.2), labels = NA, lwd = 0, lwd.ticks = 1)
# for (i in 1:length(angles)) {
#   start <- sum(angles[1:(i-1)])
#   end <- start + angles[i]
#   theta <- seq(start, end, length.out = 100) / 180 * pi
#   x <- c(0, radius[i] * cos(theta), 0)
#   y <- c(0, radius[i] * sin(theta), 0)
#   polygon(x, y, col = i, border = "white")
# }
# legend("topright", legend = labels, fill = 1:length(labels), bty = "n")
esetSel
cl3$cluster
assayMatGroupZ
plot(0,0,type='n',xlim=c(1,4),ylim=c(-1,1),axes=F,xlab="",ylab="")
clusterMeanMat<-data.frame(matrix(NA,4,4))
colnames(clusterMeanMat)<-colnames(assayMatGroupZ)
for(j in 1:4){
  assayGroupSel<-assayMatGroupZ[,j]
  for(i in 1:4){
    clusterMeanMat[i,j]<-mean(assayGroupSel[names(cl3$cluster[cl3$cluster==i])])
  }
}
# lines(seq(1,4),clusterMeanMat[1,])
# lines(seq(1,4),clusterMeanMat[2,])
# lines(seq(1,4),clusterMeanMat[3,])
# lines(seq(1,4),clusterMeanMat[4,])
colPie<-RColorBrewer::brewer.pal(4,"Pastel2")
pdf("Cluster4mfuzz_pie.pdf",width=16,height = 4)
par(mfrow=c(1,4),mar=c(3,3,3,3))
for(j in 1:4){
  plot(0,0,type='n',xlim=c(-1,1),ylim=c(-1,1),axes=F,xlab="",ylab="",main=colnames(clusterMeanMat)[j])
  polygon(sqrt(0.5)*cos(seq(0, 360, length.out = 100) / 180 * pi), sqrt(0.5)* sin(seq(0, 360, length.out = 100) / 180 * pi),col = NA, lty = 3)
  polygon(1*cos(seq(0, 360, length.out = 100) / 180 * pi),1 * sin(seq(0, 360, length.out = 100) / 180 * pi),col = NA, lty = 3)
  for(indexIn in seq(1,4)){
    angles <- seq(0,-360,length.out=indexMax+1)
    start = angles[indexIn]-15
    end = angles[indexIn+1]+15
    radius <-sqrt((clusterMeanMat[indexIn,j]+valueMax)/(2*valueMax))
    theta <- seq(start, end, length.out = 100) / 180 * pi
    x <- c(0, radius * cos(theta), 0)
    y <- c(0, radius * sin(theta), 0)
    polygon(x, y,col=colPie[indexIn],border=NA)
    abline(v = 0, col = "gray", lty = "dashed")
    abline(h = 0, col = "gray", lty = "dashed")
  }
}
dev.off()


plot(0,0,type='n',xlim=c(-1,1),ylim=c(-1,1),axes=F,xlab="",ylab="")
# plot.window(c(-1, 1), c(-1, 1))
# # axis(1, at = seq(-1, 1, by = 0.2), labels = NA, lwd = 0, lwd.ticks = 1)
# # axis(2, at = seq(-1, 1, by = 0.2), labels = NA, lwd = 0, lwd.ticks = 1)
# for (i in 1:length(angles)) {
#   start <- sum(angles[1:(i-1)])
#   end <- start + angles[i]
#   theta <- seq(start, end, length.out = 100) / 180 * pi
#   x <- c(0, radius[i] * cos(theta), 0)
#   y <- c(0, radius[i] * sin(theta), 0)
#   polygon(x, y, col = i, border = "white")
# }
# legend("topright", legend = labels, fill = 1:length(labels), bty = "n")
polyFan<-function(valueIn,indexIn,valueMax,indexMax){
  # valueIn<- -1.3895237
  # indexIn<-1
  # valueMax<-1.4
  # indexMax<-4
  plot(0,0,type='n',xlim=c(-1,1),ylim=c(-1,1),axes=F,xlab="",ylab="")
  for(indexIn in seq(1,4)){
    angles <- seq(0,-360,length.out=indexMax+1)
    start = angles[indexIn]
    end = angles[indexIn+1]
    radius <- sqrt((clusterMeanMat[indexIn,1]+valueMax)/(2*valueMax))
    theta <- seq(start, end, length.out = 100) / 180 * pi
    x <- c(0, radius * cos(theta), 0)
    y <- c(0, radius * sin(theta), 0)
    polygon(x, y,col=indexIn)
  }
}

# ● Choose a KEGG pathway topology and visualize alongside the major pathways #### 
library(xml2)
doc<-read_xml("hsa05210.xml")
root_node<-xml_root(doc)
node_list <- xml_find_all(root_node, "//entry")

graph_list<-xml_find_all(root_node,"//graphics")
idName<-typeLoc<-xLoc<-yLoc<-nameLoc<-c()
for (i in 1:length(graph_list)) {
  idName[i]<-xml_attr(node_list[i], "id")
  nameLoc[i]<-xml_attr(graph_list[i], "name")
  xLoc[i]<-as.numeric(xml_attr(graph_list[i], "x"))
  yLoc[i]<-as.numeric(xml_attr(graph_list[i], "y"))
  typeLoc[i]<-xml_attr(node_list[i], "type")
  print(i)
}
names(xLoc)<-names(yLoc)<-idName
pathName<-nameLoc[typeLoc=="map"]
pathLocx<-xLoc[typeLoc=="map"]
pathLocy<-yLoc[typeLoc=="map"]

geneName<-nameLoc[typeLoc=="gene"]
geneLocx<-xLoc[typeLoc=="gene"]
geneLocy<-yLoc[typeLoc=="gene"]

relation_list<-xml_find_all(root_node,"//relation")

src<-dst<-relations<-c()
for(i in 1:length(relation_list)){
  relations[i]<-xml_attr(xml_find_all(relation_list[i],"subtype"),"name")
  src[i]<-xml_attr(relation_list[i],"entry1")
  dst[i]<-xml_attr(relation_list[i],"entry2")
}
angleAll<-c(30,30,30,30,30,90)
ltyAll<-c(1,2,3,4,5,1)
names(ltyAll)<-names(angleAll)<-unique(relations)


#mean punch abundance
# punchP1MeanX<-apply(cbind(punchLocs[[1]][,1],punchLocs[[2]][,1],punchLocs[[3]][,1]),1,mean)
# punchP1MeanY<-apply(cbind(punchLocs[[1]][,2],punchLocs[[2]][,2],punchLocs[[3]][,2]),1,mean)
punchNames<-sapply(row.names(punchLocs[[1]]),toPlotPar)[2,]

proMat_P1_MeanPunch<-proMat_P1[,1:length(punchNames)]
proMat_P1_MeanPunch[,]<-NA
colnames(proMat_P1_MeanPunch)<-paste0("P1",punchNames)
for(j in 1:length(punchNames)){
  proMat_P1_MeanPunch[,j]<-apply(proMat_P1[,grepl(punchNames[j],colnames(proMat_P1))],1,mean,na.rm=T)
}

{
  pdf("PathPlot_MeanPunch1_CRC.pdf",width=4.90,height=3.61)
  par(mar=c(0,0,0,0))
  plot(0,0,type='n',xlim=c(200,1700),ylim=c(3600,400),axes=F,xlab="",ylab="")
for(i in 1:length(geneName)){
  text(geneLocx[i]*1.5,geneLocy[i]*3.5,labels = unlist(strsplit(geneName[i],", "))[1],pos=1,cex=0.42)
}
for(i in 1:length(relations)){
  arrows(x0=xLoc[src[i]]*1.5,x1=xLoc[dst[i]]*1.5,
         y0=yLoc[src[i]]*3.5,y1=yLoc[dst[i]]*3.5
         ,length=0.1,code = 2,cex=0.42
         ,angle=angleAll[relations[i]]
         ,lty=ltyAll[relations[i]])
}
  


# GnMatch<-read.delim("../UniProtGnMatch.tsv",stringsAsFactors = F,sep='\t',row.names = 1,header = T)
# for(i in 1:nrow(GnMatch)){
# GnMatch$Gene.Names.Main[i]<-unlist(strsplit(GnMatch$Gene.Names[i]," "))[1]
# }
uniName<-c()
for(i in 1:length(geneName)){
  gnName_sp<-unlist(strsplit(geneName[i],", "))
  mainSel<-which(GnMatch$Gene.Names.Main==gnName_sp[1])[1]
  if(length(mainSel)>0){
    uniName[i]<-row.names(GnMatch)[mainSel]
  }
  else{
    uniName[i]<-NA
  }
}
uniNameSel<-uniName[uniName%in%row.names(proMat)]
geneLocySel<-geneLocy[uniName%in%row.names(proMat)]
geneLocxSel<-geneLocx[uniName%in%row.names(proMat)]

#proMat_P1


for (i in 1:length(uniNameSel)){
  # proSel=proMat_P1#proMatList[[1]]
  proSel=proMat_P1_MeanPunch
  proSelZ<-scale(as.numeric(proSel[uniNameSel[i],]))[,1]
  # imgHeight=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[1],"_"))[2],"-"))[2])
  # imgWidth=as.numeric(unlist(strsplit(unlist(strsplit(slideContourFiles[1],"_"))[2],"-"))[1])
  # plot(0,0,type='n',xlim=c(0,imgWidth),ylim=c(imgHeight,0),axes=F,xlab="",ylab=""
  #      ,main="")#row.names(proMat)[m]
  polygon(geneLocxSel[i]*1.5+contourLines[[1]][,1]/30-40,
          geneLocySel[i]*3.5+contourLines[[1]][,2]/7+140,
          col="#9195BA",border=NA)
  pchAll<-rep(21,length(proSelZ))
  pchAll[is.na(proSelZ)]<-21
  colAll<-colGenMax1(proSelZ,heatmap_colors1)
  colAll[is.na(colAll)]<-NA
  points(geneLocxSel[i]*1.5+punchLocs[[1]][,1]/30-40,
         geneLocySel[i]*3.5+punchLocs[[1]][,2]/7+140,
         cex=0.25,pch=pchAll,bg=colAll,col=NA)
}
dev.off()
}

#   gnName_sp<-unlist(strsplit(geneName[i],", "))
#   for(m in 1:length(gnName_sp)){
#     GnMatch$Gene.Names[which(grepl(gnName_sp[1],GnMatch$Gene.Names))]
#   }
# }

#IPA####
uniNameSel<-read.csv("CCPC_IPAUniGn.csv",stringsAsFactors = F,header = F)[,2]
geneName<-read.csv("CCPC_IPAUniGn.csv",stringsAsFactors = F,header = F)[,1]

{
  pdf("PathPlot_dark_IPA.pdf",width=4.30,height=4.5)
  par(mar=c(0,0,0,0))
  plot(0,0,type='n',xlim=c(0,1700),ylim=c(3600,0),axes=F,xlab="",ylab="")
  
  # uniNameSel<-uniName[uniName%in%row.names(proMat)]
  # geneLocySel<-geneLocy[uniName%in%row.names(proMat)]
  # geneLocxSel<-geneLocx[uniName%in%row.names(proMat)]
  
  for(i in 1:length(geneName)){
    text(i%%5 *350, i%/%5 *560,labels = geneName[i],pos=1,cex=0.42)
  }
  
  for (i in 1:length(uniNameSel)){
    proSel=proMatList[[1]]
    proSelZ<-scale(as.numeric(proSel[uniNameSel[i],]))[,1]
    # polygon(geneLocxSel[i]*1.5+contourLines[[1]][,1]/30-40,
    #         geneLocySel[i]*3.5+contourLines[[1]][,2]/7+140,
    #         col="#684784",border=NA)
    polygon(i%%5 *350+contourLines[[1]][,1]/10-40,
            i%/%5 *560+contourLines[[1]][,2]/4+180,
            col="#9195BA",border=NA)
    
    pchAll<-rep(21,length(proSelZ))
    pchAll[is.na(proSelZ)]<-21
    colAll<-colGenMax(proSelZ,heatmap_colors)
    colAll[is.na(colAll)]<-NA
    points(i%%5 *350+punchLocs[[1]][,1]/10-40,
           i%/%5 *560+punchLocs[[1]][,2]/4+180,
           cex=0.22,pch=pchAll,bg=colAll,col=NA)
  }
  dev.off()
}


#volcano plots####
proMat_P1_CC
proMat_P1_PC

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
fc_P1<-fc
p_P1<-p
plot(log2(fc_P1),-log10(p_P1),xlim=c(-2,2))
abline(h=-log10(0.05),lty=2)
abline(v=log2(2^0.5),lty=2)
abline(v=log2(2^-0.5),lty=2)

p_P1_adjust<-p.adjust(p_P1,"BH")

up_P1<-row.names(proMat_P1_CC)[which(fc_P1>2^0.5&p_P1<0.05)]#31
dw_P1<-row.names(proMat_P1_CC)[which(fc_P1<2^-0.5&p_P1<0.05)]#4

names(fc_P1)<-names(p_P1)<-row.names(proMat_P1_CC)

pdf("Volco_CCPC.pdf")
plot(log2(fc_P1),-log10(p_P1),xlim=c(-2,2),pch=19,col="grey",xlab="log2 (Fold change)"
     ,ylab="log10 p value")
GnMatch[up_P1,3]
points(log2(fc_P1[up_P1]),-log10(p_P1[up_P1]),pch=19,cex=1.6,col="#4B0055")
points(log2(fc_P1[dw_P1]),-log10(p_P1[dw_P1]),pch=19,cex=1.6,col="#274983")
text(log2(fc_P1[up_P1]),-log10(p_P1[up_P1])+0.25,GnMatch[up_P1,3])
text(log2(fc_P1[dw_P1]),-log10(p_P1[dw_P1])+0.25,GnMatch[dw_P1,3])
dev.off()
# abline(h=-log10(0.05),lty=2)
# abline(v=log2(2^0.25),lty=2)
# abline(v=log2(2^-0.25),lty=2)
outBind<-cbind(fc,p)

updw_P1_ind<-c(which(fc_P1>2^0.3&p_P1<0.05),which(fc_P1<2^-0.3&p_P1<0.05))


row.names(outBind)<-row.names(proMat_P1_CC)
write.csv(outBind[updw_P1_ind,],"CCPC_fcp_0.3.csv")

write.csv(up_P1,"CCup_P1.csv")
write.csv(dw_P1,"CCdw_P1.csv")
# proMat_PC<-proMat[,grepl("PC",colnames(proMat))]
# proMat_CC<-proMat[,grepl("CC",colnames(proMat))]
# 
# fc<-p<-c()
# for(i in 1:nrow(proMat_PC)){
#   sel1<-as.numeric(proMat_PC[i,])
#   sel2<-as.numeric(proMat_CC[i,])
#   if(sum(!is.na(sel1))>2&sum(!is.na(sel2))>2){
#     fc[i]<-mean(sel1,na.rm=T)/mean(sel2,na.rm=T)
#     p[i]<-t.test(sel1,sel2)$p.value
#   }else{
#     fc[i]<- p[i]<-NA
#   }
# }
# plot(log2(fc),-log10(p),xlim=c(-2,2))
# abline(h=-log10(0.05),lty=2)
# abline(v=log2(2^0.5),lty=2)
# abline(v=log2(2^-0.5),lty=2)
# 
# up<-row.names(proMat_PC)[which(fc>2^0.5&p<0.05)] #26
# dw<-row.names(proMat_PC)[which(fc<2^-0.5&p<0.05)] #40

#sankey plots####
SankIn<-list.files("Sankey/")
path_pro_cl<-c()
for(i in 1:length(SankIn)){
  readIn<-read.delim(paste0("Sankey/",SankIn[i]),stringsAsFactors = F,sep='\t',header = T,fill = T)
  readInSel<-readIn[readIn$X.log.p.value.> -log10(0.02),]
  for(m in 1:nrow(readInSel)){
    pro_sp<-unlist(strsplit(readInSel$Molecules[m],","))
    for(n in 1:length(pro_sp)){
      path_pro_cl<-c(path_pro_cl,paste(readInSel$Ingenuity.Canonical.Pathways[m],pro_sp[n],i,sep=';'))
    }
  }
}

df_proPath<-data.frame(matrix(NA,nrow=length(path_pro_cl),ncol=3),stringsAsFactors = F)
df_proPath[,1]<-sapply(sapply(path_pro_cl,strsplit,";"),"[[",2)
df_proPath[,2]<-sapply(sapply(path_pro_cl,strsplit,";"),"[[",1)
df_proPath[,3]<-sapply(sapply(path_pro_cl,strsplit,";"),"[[",3)

df_proPath[df_proPath[,3]==4,3]<-0
df_proPathOd<-df_proPath[order(df_proPath[,3]),]

library(ggplot2)
library(ggalluvial)
data <- data.frame(
  from = (df_proPathOd[,1]),
  to = (df_proPathOd[,2]),
  value = rep(1,nrow(df_proPathOd))
)
data$from<-factor(data$from, levels = unique(data$from))
data$to<-factor(data$to, levels = unique(data$to))


pdf("Sankey.pdf")
ggplot(data, aes(axis1 = from, axis2 = to, y = value)) +
  geom_alluvium(aes(fill = from), width = 1/12) +
  geom_stratum(width = 1/12, fill = "gray", color = "black") +
  theme_minimal()
dev.off()

pdf("Sankey_legend.pdf",width=23)
par(mfrow=c(1,2))
b1<-barplot(as.matrix(table(data$from)),beside = F,col=0)
text(b1,cumsum(table(data$from))-0.5,names(table(data$from)))

b2<-barplot(as.matrix(table(data$to)),beside = F,col=0)
text(b2,cumsum(table(data$to))-0.5,names(table(data$to)))
dev.off()

ggplot(data, aes(axis1 = from, axis2 = to, y = value, fill = to)) +
  geom_alluvium(aes(fill = to), width = 1/12) +
  geom_stratum(width = 1/12, color = "black") +
  scale_fill_manual(values = c("D" = "red", "E" = "blue", "F" = "green", "G" = "orange"),
                    labels = c("Label D", "Label E", "Label F", "Label G")) +
  theme_minimal()
