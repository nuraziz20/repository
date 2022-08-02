####IMPORTANT! Please download the github version to use the most recent version of TCGAbiolinks
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
###Use devtools::install_github(repo = "ELELAB/TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(recount)
###conversion of uuids to TCGA barcodes
library(TCGAutils)
library(limma)
library(biomaRt)


####Please download the github version to use the right version of TCGAbiolinks and TCGAutils

###Use devtools::install_github()


####Function to convert ENSG to Symbol
convert.ENSG.Symbol<-function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  #merge(DEG.ucs,G_list,by.x="gene",by.y="ensembl_gene_id")
  return(G_list)
}



########Query from Recount2 platform: SKC Carcinoma#######
skcm.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue="skin")
skcm.recount.tcga<-TCGAquery_recount2(project="TCGA", tissue="skin")

#to get the SE object
SE.skcm.recount.gtex <- skcm.recount.gtex$GTEX_skin
SE.skcm.recount.tcga <- skcm.recount.tcga$TCGA_skin

# to get the read-count matrix only
#matrix <- assays(SE.ucs.recount.gtex)$counts


#### The steps below are needed to have the right correspondance beetween barcodes (TCGA) and UUID (recount)
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")


query.skcm<- GDCquery(project = "TCGA-SKCM",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts")


samplesDown.skcm <- getResults(query.skcm,cols=c("cases"))


###tumor samples for skin cancer
dataSmTP.skcm <- TCGAquery_SampleTypes(barcode = samplesDown.skcm,
                                      typesample = "TP")

###to check that there are no NT samples
dataSmNT.skcm <- TCGAquery_SampleTypes(barcode = samplesDown.skcm,
                                      typesample = "NT")


#####Preparing/scaling Recount2 data because it was sequenced using Rail-RNA
BiocManager::install(c("recount"))
library(recount)

eset.gtex<-assays(scale_counts(skcm.recount.gtex$GTEX_skin, round = TRUE))$counts
eset.tcga<-assays(scale_counts(skcm.recount.tcga$TCGA_skin, round = TRUE))$counts

#### Check that the number of reads is less than or equal to 40 million
rse_scaled <- scale_counts(skcm.recount.gtex$GTEX_skin, round = TRUE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6

###replacing UUIDs with TCGA barcodes:
colnames(eset.tcga)<-colData(skcm.recount.tcga$TCGA_skin)$gdc_cases.samples.portions.analytes.aliquots.submitter_id

###Removing version (number after ".")
rownames(eset.gtex) <- gsub("\\..*", "", rownames(eset.gtex))
rownames(eset.tcga) <- gsub("\\..*", "", rownames(eset.tcga))

####Segregate between primary tumors and normal samples
eset.tcga.cancer<-eset.tcga[,which(colData(skcm.recount.tcga$TCGA_skin)$gdc_cases.samples.sample_type=="Primary Tumor")]
eset.tcga.normal<-eset.tcga[,which(colData(skcm.recount.tcga$TCGA_skin)$gdc_cases.samples.sample_type=="Solid Tissue Normal")]
####


##merging data by row names
dataPrep.skcm<-merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by=0, all=TRUE)

rownames(dataPrep.skcm)<-dataPrep.skcm$Row.names
dataPrep.skcm$Row.names<-NULL


dataNorm.skcm <- TCGAanalyze_Normalization(tabDF = dataPrep.skcm,
                                          geneInfo = geneInfoHT,
                                          method = "gcContent")

# Error in `.rowNamesDF<-`(x, value = value) : 
#  duplicate 'row.names' are not allowed

## delete duplicated row
dataFilt.skcm <- TCGAanalyze_Filtering(tabDF = distinct(dataPrep.skcm),
                                      method = "quantile", 
                                      qnt.cut =  0.25)
colnames(dataFilt.skcm)

##converting ensenmbl gene ids to huugo gymbols using Biomart package
conversion.table<-convert.ENSG.Symbol(rownames(dataFilt.skcm))

conversion.inter<-intersect(conversion.table[-which(conversion.table$hgnc_symbol==""),]$ensembl_gene_id, rownames(dataFilt.skcm))
conversion.table2<-conversion.table[which(conversion.table$ensembl_gene_id %in% conversion.inter),]

#non-unique values when setting 'row.names': ‘ENSG00000230417’, ‘ENSG00000276085’
##These genes were removed because they cause duplicated rows
c("ENSG00000230417", "ENSG00000276085")
dataFilt.skcm<-dataFilt.skcm[-which(rownames(dataFilt.skcm) %in% c("ENSG00000230417", "ENSG00000276085")), ] #run again for conversion table
rownames(conversion.table2)<-conversion.table2$ensembl_gene_id

conversion.table[-which(conversion.table$hgnc_symbol==""),]
dataFilt.skcm.hugo<-dataFilt.skcm[conversion.inter,]

##after removing duplicate row names, we merge the hugo symboles with DEGs table
dataFilt.skcm.hugo<-merge(dataFilt.skcm.hugo, conversion.table2, by=0)
rownames(dataFilt.skcm.hugo)<-dataFilt.skcm.hugo$Row.names
dataFilt.skcm.hugo$Row.names<-NULL


write.csv(dataFilt.skcm.hugo,"GTEx TCGA tumor skcm.csv")
dataFilt.skcm.hugo
selected<-dataFilt.skcm.hugo %>% 
  filter(hgnc_symbol=="LAMC2")

t.selected<-as.data.frame(t(selected))
colnames(t.selected)[1]<-"exp"
t.selected<-as.data.frame(t.selected) %>% 
  mutate(Group=substr(rownames(t.selected),1,1))
head(t.selected)

# Plotting
colnames(t.selected)
#remove e(ensemble) and h(hgsc) group
t.selected<-t.selected %>% 
  filter(Group == 'S' | Group =='T') %>% 
  mutate(Group=case_when(Group == "S" ~ "GTEx Normal",
                      Group == "T" ~ "TCGA-SKCM Tumor"))

ggplot(data=t.selected,aes(x=Group,y=log2(as.numeric(`exp`)),fill=Group))+
  stat_boxplot(geom ='errorbar',width=.2) 
  geom_jitter(colour='gray90',alpha=.8,width = .2)+
  geom_boxplot(width=.6)+
  theme_Publication()+
  scale_fill_aziz()

# statistics
wilcox.test(log2(as.numeric(`exp`)) ~ Group, data = t.selected, exact = FALSE)
t.selected %>% 
  count(Group)


#### Many genes #####
selected<-dataFilt.skcm.hugo %>% 
  filter(hgnc_symbol=="METTL18"|
         hgnc_symbol=="CAMKMT"|
         hgnc_symbol=="EEF2KMT"|
         hgnc_symbol=="ETFBKMT"|
         hgnc_symbol=="METTL21A"|
         hgnc_symbol=="EEF1AKMT3"|
         hgnc_symbol=="METTL21C"|
         hgnc_symbol=="VCPKMT"|
         hgnc_symbol=="METTL22"|
         hgnc_symbol=="METTL23")

rownames(selected)<-selected$hgnc_symbol

t.selected<-as.data.frame(t(selected))
colnames(t.selected)
t.selected<-as.data.frame(t.selected) %>% 
  mutate(Group=substr(rownames(t.selected),1,1))
head(t.selected)

# Plotting
colnames(t.selected)
#remove e(ensemble) and h(hgsc) group
t.selected<-t.selected %>% 
  filter(Group == 'S' | Group =='T')

#apply log2 to 10 column
m16<-t.selected %>% 
  mutate(log2.METTL18 = log2(as.numeric(METTL18)),
         log2.METTL23 = log2(as.numeric(METTL23)),
         log2.METTL21A = log2(as.numeric(METTL21A)),
         log2.CAMKMT = log2(as.numeric(CAMKMT)),
         log2.METTL21C = log2(as.numeric(METTL21C)),
         log2.ETFBKMT = log2(as.numeric(ETFBKMT)),
         log2.EEF1AKMT3 = log2(as.numeric(EEF1AKMT3)),
         log2.EEF2KMT = log2(as.numeric(EEF2KMT)),
         log2.VCPKMT = log2(as.numeric(VCPKMT)),
         log2.METTL22 = log2(as.numeric(METTL22))) %>% 
  dplyr::select(Group, 
         log2.METTL18 ,
         log2.METTL23,
         log2.METTL21A,
         log2.CAMKMT,
         log2.METTL21C,
         log2.ETFBKMT,
         log2.EEF1AKMT3,
         log2.EEF2KMT,
         log2.VCPKMT,
         log2.METTL22)

m16<-m16 %>% 
  rename(METTL18=log2.METTL18 ,
         METTL23=log2.METTL23,
         METTL21A=log2.METTL21A,
         CAMKMT=log2.CAMKMT,
         ETFBKMT=log2.ETFBKMT,
         EEF1AKMT3=log2.EEF1AKMT3,
         EEF2KMT=log2.EEF2KMT,
         VCPKMT=log2.VCPKMT,
         METTL22=log2.METTL22) %>% 
  gather('Gene','Value',c(METTL18,
                          METTL23,
                          METTL21A,
                          CAMKMT,
                          ETFBKMT,
                          EEF1AKMT3,
                          EEF2KMT,
                          VCPKMT,
                          METTL22)) %>% 
  dplyr::select(Gene,Value,Group) 


ggplot(data=m16,aes(x=Group,y=Value,fill=Group))+
  stat_boxplot(geom ='errorbar',width=.2) +
  geom_jitter(colour='gray90',alpha=.8,width = .2)+
  geom_boxplot(width=.6)+
  facet_grid(~Gene)+
  theme_Publication()+
  scale_fill_aziz()


ggsave('skin.png')

head(m16)
# statistics
kruskal.test(Value ~ , data = m16, exact = FALSE)

m16 %>% 
  dplyr::filter(Gene=="EEF2KMT") %>% 
  kruskal.test(EEF2KMT ~ Group, data = m16)

t.selected %>% 
  count(Group)
