

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggplot2")
library(ggplot2)

##读入数据#########################
#biowolf<-read.table(file="DEGs_information.txt",header=TRUE,sep="\t",row.names = 1)
biowolf<-read.csv(file="DEGs_information.csv",header=TRUE,sep=",",row.names = 1)
colnames(biowolf)
biowolf$color <- ifelse(biowolf$P.Value<0.05 & abs(biowolf$logFC)>= 0,ifelse(biowolf$logFC> 0,'red','blue'),'gray')
#??logFC??ֵ
color <- c(red = "red",gray = "gray",blue = "blue")
pdf(file="DEGs_MA_P0.05.pdf",width=10,height=8,onefile=FALSE)
p <- ggplot(biowolf, aes(logFC, -log10(P.Value), col = color))+
  xlim(-3,3)+ylim(0,6)+geom_point()+theme_bw()+scale_color_manual(values = color)+
#??DEGs?ĺ?????????????ֵ????Сֵ??��??xlim??ylim??��??֤???̵???ͼ????ʾ  
  labs(title="",x="log2 (fold change)",y="-log10 (P.Value)")+
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6)+
  geom_vline(xintercept = 0, lty=4,col="grey",lwd=0.6)+
  theme(plot.title = element_text(hjust = 0.5),
        #legend.position = "none",
        #legend.position = "topright",
        #legend.position =  c(0.95, 0.9),
        #legend=c("Up","Down","nonsignificant"),
        #col=c("red","blue","gray"),
        #legend.margin= T,
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#xinterceptֵ??logFC?????õ?ֵ??Ӧ??
print(p)
dev.off()


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")  
BiocManager::install("vsn") 
BiocManager::install("pheatmap") 

library("vsn")
library("pheatmap")

##?????????ȱ???��????????
data_FPKM <- read.table("sampleExp.txt",sep="\t",header=T,check.names = F)   ###用的是基因表达量矩阵，，，横为样本，纵为基因???
rownames(data_FPKM) <- data_FPKM[,1]
data_FPKM <- data_FPKM[,-1]
#data_FPKM<-log2(data_FPKM)

##?????????????????Ľ???
input_data<-"DEGs_information.csv"   ##用的数据是limma差异表达分析跑出来的结果，包括gene，Pvalue,adj p等等
data <-read.table(file=input_data,header=T,sep=",")
p <-0.05
logFoldChange <- 0    ###看情况改

res_de <- subset(data,data$P.Value < p)

res_de_up <- subset(res_de, res_de$logFC>=logFoldChange)
res_de_up_top50_id <- as.vector(head(res_de_up$X,50))
res_de_dw <- subset(res_de, res_de$logFC<=-logFoldChange)
res_de_dw_top50_id <- as.vector(head(res_de_dw$X,50))

red_de_top100 <- c(res_de_up_top50_id, res_de_dw_top50_id)
#red_de_top100 <- c(res_de_up, res_de_dw)
#当上下调差异表达基因不够1001个则用上一行进行跑。

##data_FPKM <- count(data_FPKM, normalized=TRUE)
##data_FPKM_map<-apply(data_FPKM, 1, mad)
##data_FPKM <- data_FPKM[order(data_FPKM_map, decreasing=T), ]
red_de_top100_expr <- data_FPKM[rownames(data_FPKM) %in% res_de_dw_top50_id,]


Type<-"group.txt"   ##用的数据是limma差异表达分析跑出来的结果，包括gene，Pvalue,adj p等等
ann <-read.table(file=Type,header=T,sep="\t",check.names = F)
group <- t(ann)
group <- data.frame(group)
#library(tidyverse)
#select(ann1,group=1)

red_de_top100_expr <- red_de_top100_expr[,rownames(group)]

identical(rownames(group),colnames(red_de_top100_expr))
pdf(file="DEGs_pheatmap_down_50.pdf",width=9,height=6,onefile=FALSE)

#annotation_col=data.frame(group=rep(c("test","control"),c(6,6)))
#需要改样本标签和对应样本数。
#rownames(annotation_col)=colnames(data_FPKM)

#因为red_de_top100_expr样本名是用-连接，
#所以需要把表达量数据（data_FPKM）和分组文件（annotation_col）中样本名的.替换成-
#colnames(data_FPKM)<-gsub("\\.","\\-",colnames(data_FPKM))


##2.##样本顺序不一样时,前期应改变样本顺序就可避免这类问题

pheatmap(red_de_top100_expr,
         annotation_col = group,
         border_color=NA,
         cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row=7,show_colnames=F,scale="row",
         cluster_cols=FALSE,main="down DEGs heatmaps")
dev.off()




