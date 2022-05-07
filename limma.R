
library(edgeR)
library(limma)

data_FPKM <- read.table("sampleExp.txt",sep="\t",header=T,check.names = F)
rownames(data_FPKM) <- data_FPKM[,1]
data_FPKM <- data_FPKM[,-1]

clinical<-read.table("group.txt",sep="\t",header=F)
clinical<-t(clinical)
clinical<-data.frame(clinical)
colnames(clinical)<-c("var1","var2")
rownames(clinical)<-clinical$var1

clinical<-clinical[as.vector(colnames(data_FPKM)),] 
group<-factor(clinical$var2)
design <- model.matrix(~0+group, data=clinical)
colnames(design) <- c(levels(group))
identical(colnames(data_FPKM),rownames(clinical))

fit<-lmFit(data_FPKM,design)
contrast.matrix <- makeContrasts("control-case",levels=design) 
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
result1<-topTable(fit2,coef=1,adjust="BH",number=dim(data_FPKM)[1])
write.csv(result1,file="DEGs_information.csv",quote=F,row.names=T)

logFc=0
PValue=0.05
diffSig = result1[(result1$P.Value < PValue & (result1$logFC>logFc | result1$logFC<(-logFc))),]
diffSig <- cbind(gene=rownames(diffSig),diffSig)
write.table(diffSig, file="diffSig_p0.05.xls",row.names=F,sep="\t",quote=F)

#UP
diffUp = result1[(result1$P.Value < PValue & (result1$logFC>logFc)),]
diffUp <- cbind(gene=rownames(diffUp),diffUp)
write.table(diffUp, file="up_p0.05.xls",row.names=F,sep="\t",quote=F)

#Down
diffDown = result1[(result1$P.Value < PValue & (result1$logFC<(-logFc))),]
diffDown <- cbind(gene=rownames(diffDown),diffDown)
write.table(diffDown, file="down_p0.05.xls",row.names=F,sep="\t",quote=F)
