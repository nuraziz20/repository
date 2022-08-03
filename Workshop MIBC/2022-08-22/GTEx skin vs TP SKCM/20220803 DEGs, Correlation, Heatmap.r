# library
install.packages("tidyverse") #omit this if you have install tidyverse/dplyr
install.packages("dplyr")

library(tidyverse)
library(dplyr)

#### import dataset SKCM , and name it: skcm (remember set heading to "yes")
skcm<-GTEx.TCGA.tumor.skcm

# check dimension
dim(skcm) #it should show : 30831 (row of genes) and 1080 (column name of sample/tissue of patients)

View(skcm) # you notice that there is colnames ensemble (column 1 and 1079) and hgnc (column 1080)

#lets move to the fisrt and second column

skcm<-skcm %>% 
  select(hgnc_symbol,2:1078) #this will select column hgnc_symbol and samples

# now we will remove the normal sample
tumor.skcm<-skcm %>% 
  select(hgnc_symbol,contains("TCGA"))

# prepare your deg (you can get from rnaseq from ebiogen/public rnaseq dataset)
# import example DEGs and name it : DEGs
DEGs<-DEGs.example

# lets merge two dataframe by left join
# first make sure the most left column name is same
colnames(tumor.skcm)[1]
colnames(DEGs)
#renaming
colnames(tumor.skcm)[1]<-"Gene"

DEGs<-DEGs %>% 
  left_join(tumor.skcm)
View(DEGs)

#remove ICT1 gene (not detected using this name)
filtered.DEGs<-DEGs %>% 
  drop_na()

#### gene-gene correlation ####
# flip the table
t.DEGs<-t(filtered.DEGs)
colnames(t.DEGs)<-t.DEGs[1,]
t.DEGs<-t.DEGs[-1,]
char.DEGs<-t.DEGs
m.DEGs<-matrix(as.numeric(char.DEGs),
               ncol=ncol(char.DEGs))
df.DEGs<-as.data.frame(m.DEGs)
colnames(df.DEGs)<-colnames(char.DEGs)
rownames(df.DEGs)<-rownames(char.DEGs)

#convert log
df.DEGs<-log(df.DEGs)

# 2 gene
cor.test(df.DEGs$METTL23,df.DEGs$MRPS7)
df.DEGs %>% 
  ggplot(aes(x=METTL23,y=MRPS7))+
  geom_smooth(method=lm, fullrange=FALSE,se=TRUE)+
  geom_point()

# make table correlation value
library(correlation)

correlation::correlation(df.DEGs,
                         include_factors = TRUE, method = "auto") %>% 
  filter(Parameter1 == 'METTL23')

# make graph correlation METTL23 with other genes in one plot
df.DEGs %>% 
  gather('Gene','Value',2:9) %>% 
  select(Gene,Value,METTL23) %>% 
  ggplot(aes(x=Value,y=METTL23))+
  geom_smooth(method=lm, fullrange=FALSE,se=TRUE)+
  geom_point()+
  facet_grid(~Gene,scales = "free_x")+
  theme(aspect.ratio = 1)

# matrix correlation heatmap
library(ggcorrplot)
library(ggstatsplot)
ggstatsplot::ggcorrmat(
  data = df.DEGs,
  type = "parametric", # parametric for Pearson, nonparametric for Spearman's correlation
  colors = c("steelblue","white","darkred"), # change default colors
  title = "Correlation matrix of DEGs",
  subtitle = "Skin Cutaneous Melanoma (SKCM-TCGA)",
  matrix.type ='upper',
  ggcorrplot.args = list(outline.color = "white", 
                         hc.order = TRUE, #clustering
                         pch.cex=3, # x size
                         lab_size=3.25) #label size
  
)+ggplot2::theme(aspect.ratio = 1,
                 axis.text = element_text(size=10, colour = 'black',family = 'Arial'),
                 axis.text.x= element_text(family = 'Arial',hjust =0,vjust = 1)
)+scale_x_discrete(position = 'top')

# heatmap each patients
library(pheatmap)
library(viridis)
library(ComplexHeatmap)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

hm.DEGs <- apply(df.DEGs, 2, cal_z_score)
hm.DEGs[hm.DEGs < -2] = -2
hm.DEGs[hm.DEGs > 2] = 2

colmy<-colorRampPalette(c("#4DBBD5B2","white","#E64B35B2"))(50)
colmy<-paste0(colmy,'7f')
ComplexHeatmap::pheatmap(t(hm.DEGs),
                         show_colnames = F,
                         color = colmy,
                         #row_title = "Genes", row_title_rot = 0,
                         column_title = "Patients",
                         heatmap_legend_param = list(title = gt_render("<span style='color:black'>*z-score*</span>")))
