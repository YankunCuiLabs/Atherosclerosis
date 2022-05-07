
install.packages("tidyr")
install.packages("ggplot2")
install.packages("tidyverse")
library(tidyverse)
library(tidyr)
library(ggplot2)
require(ggpubr)
require(ggsci)
require(cowplot)
library(RColorBrewer)

mycol <- brewer.pal(2,'Dark2')

hlaexp<-read.table("24_gene_exp.txt",sep="\t",header=T,check.names = F)
#hlaexp <- hlaexp[,-3]
rownames(hlaexp) <- hlaexp[,1]
#rownames(hlaexp) <- hlaexp$ID 
hlaexp=hlaexp[,-1]

boxdata <- gather(hlaexp,HLA,exp,2:ncol(hlaexp))
boxdata[,'exp'] <- sapply(boxdata[,'exp'],function(x){log2(x+1)})

boxdata$group <- factor(boxdata$group,levels = c('Control','Case'))


pdf('24_exp_new_no_point.pdf',width =8,height = 5)
ggboxplot(boxdata,x='HLA',y='exp',color ='group',
          ylab = 'Relative Expression',
          xlab ='',palette =mycol,
          #palette = "d3",  #nejm #jco #jama#npg#d3#aaas
          #add ='jitter',add.params = list(size=0.1,shape=16),   #add ='jitter'
          width=0.5)+
  rotate_x_text(45)+
  #coord_flip()+   #颠倒X轴和Y轴
  theme(axis.text=element_text(size=8),
        axis.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  stat_compare_means(size = 4,aes(group=group),
                     label = 'p.signif',method = 'wilcox.test',hide.ns = F)  ##检验方法可以改
dev.off()


