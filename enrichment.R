library(ggplot2)
library(stringr)

slimgo = read.table("GO_input.txt",header=T,sep="\t")
slimgo$Term=factor(slimgo$Term,levels=slimgo$Term)

pdf(file="GO.pdf",width=12,height=6,onefile=FALSE)
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5") ##自定义颜色


p<-ggplot(slimgo, aes(x=Term,y=-log10(P_Value),fill=Category))+  
  geom_bar(stat="identity", width=0.6)+
  coord_flip()+    #guides(fill=T)+         
  scale_x_discrete(limits=rev(levels(slimgo$Term)))+
  #scale_fill_manual(values = CPCOLS)+ 
  labs(x="GO Terms",y="-log10 (P_Value)",title="")+
  #scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=30))+
  geom_text(mapping = aes(label = Count),size=3,vjust=0.5 ,hjust=-0.25) 

mytheme <- theme_bw() + theme(
         axis.text=element_text(size=10,colour="black"),
        axis.title=element_text(size=10,colour="black"))
p=p+mytheme
p
dev.off()

