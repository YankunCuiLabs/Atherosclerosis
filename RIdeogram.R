
install.packages("RIdeogram")
library(RIdeogram)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


data(human_karyotype, package="RIdeogram")
#karyotype <- read.csv("karyotype.csv")
human_karyotype


data(gene_density, package="RIdeogram")
head(gene_density)

Random_RNAs_500 <- read.csv("xiaoyufu0.2_input_mark.csv")
head(Random_RNAs_500)


ideogram(karyotype = human_karyotype, 
         
         overlaid = gene_density, 
         colorset1 = c("#2c7fb8", "white", "#e34a33"),
         label = Random_RNAs_500, 
         label_type = "marker",
         width = 200, 

         Lx = 160, 
         Ly = 20, 
         output = "chromosome.svg") 

svg2pdf("chromosome.svg", 
        width = 12, height = 8, 
        dpi = 300)


sessionInfo()


