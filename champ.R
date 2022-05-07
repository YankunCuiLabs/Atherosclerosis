
library(data.table)
BiocManager::install("ChAMP",ask = F,update = F)
  
group <- read.table("group.txt",header =T,check.names = F)
b=group
rownames(b)=b[,1]   ##得到行-样本，列-分组文件

library(data.table)
a=fread("sampleExp.txt",data.table = F)  
#a=read.table("data.txt", header =T,sep="\t",check.names = F)  
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
a[1:4,1:4]
c <- na.omit(a)  ##删除整行为NA值
c[1:4,1:4]
beta=as.matrix(c)
library(impute)
beta=impute.knn(beta)
betaData=beta$data
betaData=betaData+0.00001
a=betaData


identical(colnames(a),rownames(b))
aa<-cbind(id=rownames(a),a)
write.table(aa,"treated_data.txt",sep="\t",quote=F,row.names=F)

library(ChAMP)
myLoad=champ.filter(beta = a,pd = b)
#myLoad
save(myLoad,file = 'step1-myLoad.Rdata')
 

champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Group)

#champ.QC()函数结果三张图，如有mdsPlot (Multidimensional Scaling Plot)，主要看看不同分组样本是否分开；
#densityPlot，每个样本的beta值的分布图，主要看看有无异常的样本；
#dendrogram，样本的聚类图

load(file = 'step1-myLoad.Rdata')

if(T){  
  myNorm <- champ.norm(beta=myLoad$beta,
                       rgSet=myLoad$rgSet,
                       mset=myLoad$mset,
                       resultsDir="./CHAMP_Normalization/",
                       method="BMIQ",
                       plotBMIQ=FALSE,
                       arraytype="450K",
                       cores=5)
  dim(myNorm) 
  pD=myLoad$pd
  save(myNorm,pD,file = 'step2-champ_myNorm.Rdata')
}

load(file = 'step2-champ_myNorm.Rdata')
beta.m=myNorm
group_list=myLoad$pd$Group
dim(beta.m) 


library(ChAMP)
champ.SVD(beta=myNorm,
          rgSet=NULL,
          pd=myLoad$pd,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages/")


rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-myLoad.Rdata')
load(file = 'step2-champ_myNorm.Rdata')
group_list=myLoad$pd$Group
table(group_list)


myDMP <- champ.DMP(beta = myNorm,pheno=group_list,adjPVal = 1)  
head(myDMP[[1]])   
champDiff=myDMP[[1]]  
write.csv(champDiff,"champDiff_adjp1.csv",quote=F)
save(myDMP,file = 'step3-output-myDMP.Rdata')


myDMP_adjP0.05 <- champ.DMP(beta = myNorm,pheno=group_list,adjPVal = 0.05)  
champDiff_0.05 =myDMP_adjP0.05[[1]]
write.csv(champDiff_0.05,"champDiff_adjp0.05.csv",quote=F)
