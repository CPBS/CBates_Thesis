require(stringr)
require(BiocInstaller)
require(ggplot2)
require(UpSetR)
require(dplyr)
require(ggbeeswarm)
library(patchwork)
require(forcats)
require(topGO)
biocLite("org.Sc.sgd.db")
require(org.Sc.sgd.db)
Slider = read.csv("~/Downloads/Disorder Stuff/allS288C_Prot_SliderPredScores.csv")
for(i in 1:length(Slider$Prot..ID)){
  Slider$ID[i] = substr(Slider$Prot..ID[i],3,9)
}


theme_Publication_Gridlines <- function(base_size=9, base_family="Calibri", gridlines = FALSE,legendPos="bottom",legendDir = "horizontal") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background  = element_blank(),
            plot.background = element_rect(fill="transparent", colour=NA),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = 8),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            legend.key = element_rect(fill="transparent", colour=NA),
            legend.position = legendPos,
            legend.background = element_rect(fill="transparent", colour=NA),
            legend.direction = legendDir,
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold"),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key.height=unit(1,"line")
    ))
}



#Load in all RBP datasets I could find from yeast and format them so we have just 2 columns - column 1 is the common name, column 2 is the SGD ID.
require(gdata)
RBP1=read.xls("~/Downloads/RNA_Binding_Proteins/Beckmann_2017_Hentze/ncomms10127-s2.xlsx")
RBP11=RBP1[3:length(RBP1$X),c(1,2)]
colnames(RBP11) = c("Gene","ENSEMBL")
RBP2=read.xls("~/Downloads/RNA_Binding_Proteins/Brannan_2016_Yeo/1-s2.0-S1097276516305202-mmc5_YeastHumanFly.xlsx",sheet=2)
RBP22=RBP2[,c(3,9)]
RBP22=RBP22[which(RBP22$yeast_RBP_flag==1),]
for(i in 1:length(RBP22$Gene)){
  RBP22$ENSEMBL[i]=get(as.character(RBP22[i,1]), org.Sc.sgdCOMMON2ORF)
}
RBP22=RBP22[,c(1,3)]
colnames(RBP22)=c("Gene","ENSEMBL")
RBP3=read.xls("~/Downloads/RNA_Binding_Proteins/Gonzalez_2015_Gerber/nsmb.3128-S3_UseYeastmRBPOmeColumn.xls",sheet=2)
RBP33=RBP3[2:length(RBP3$X),c(2,1,4)]
RBP33=RBP33[which(RBP33[,3]=="Y"),]
RBP33$X.1=gsub(" ","", RBP33$X.1)
RBP33=RBP33[,c(1,2)]
colnames(RBP33)=c("Gene","ENSEMBL")
RBP4=read.xls("~/Downloads/RNA_Binding_Proteins/Hogan_2008_OBrown/journal.pbio.0060255.st001.XLS",sheet=2)
RBP44=RBP4[,c(2,1)]
colnames(RBP44)=c("Gene","ENSEMBL")
RBP5=read.xls("~/Downloads/RNA_Binding_Proteins/Klass_2013_OBrown/Supplemental_Table_S4_RNADepRNABinding.xls")
RBP55=RBP5[2:length(RBP5$Systematic.Name),c(2,1)]
colnames(RBP55)=c("Gene","ENSEMBL")
RBP6=read.xls("~/Downloads/RNA_Binding_Proteins/Klass_2013_OBrown/Supplemental_Table_S5_NovelRBPs.xls")
RBP66=RBP6[2:length(RBP6$Systematic.Name),c(2,1)]
colnames(RBP66)=c("Gene","ENSEMBL")
RBP7=read.xls("~/Downloads/RNA_Binding_Proteins/Klass_2013_OBrown/Supplemental_Table_S6_KnownRBPs.xls")
RBP77=RBP7[,c(2,1)]
colnames(RBP77)=c("Gene","ENSEMBL")
#Keep only the data from Klass that is unique - remove repeated observations
KlassList = list(RBP55$Gene, RBP66$Gene, RBP77$Gene)
KlassList=setNames(KlassList,c(5,6,7))
klassMat = list_to_matrix(KlassList) #Make a comparison matrix
simpleMat = list_to_matrix(KlassList)
simpleMat =as.data.frame(simpleMat) #Convert to a dataframe
RBPKlass = data.frame(Gene=rownames(simpleMat))#Take the row names
#subset(simpleMat,`6`==1 | `5`==1 | `7`==1) #Keep those IDs which are unique
for(i in 1:length(RBPKlass[,1])){
  if(substr(as.character(RBPKlass[i,1]),nchar(as.character(RBPKlass[i,1])),nchar(as.character(RBPKlass[i,1])))=="C" | substr(as.character(RBPKlass[i,1]),nchar(as.character(RBPKlass[i,1])),nchar(as.character(RBPKlass[i,1])))=="Y" |substr(as.character(RBPKlass[i,1]),nchar(as.character(RBPKlass[i,1])),nchar(as.character(RBPKlass[i,1])))=="W"){
    RBPKlass$ENSEMBL[i]=as.character(RBPKlass[i,1])
  }
  else{
  RBPKlass$ENSEMBL[i]=get(as.character(RBPKlass[i,1]), org.Sc.sgdCOMMON2ORF)
  }
  print(i)
}
RBP8=read.xls("~/Downloads/RNA_Binding_Proteins/Scherrer_2010_Gerber/journal.pone.0015499.s002.XLS")
RBP88=RBP8[,c(2,1)]
colnames(RBP88)=c("Gene","ENSEMBL")
RBP9=read.xls("~/Downloads/RNA_Binding_Proteins/Tsvetanova_2010_OBrown/journal.pone.0012671.s001_DatabaseSearch.XLS")
RBP99=RBP9[,c(1,2)]
RBP99$Protein=gsub(" ","", RBP99$Protein)
RBP99$Gene=gsub(" ","",RBP99$Gene)
colnames(RBP99)=c("Gene","ENSEMBL")
RBP10=read.xls("~/Downloads/RNA_Binding_Proteins/Mitcehll_2013_ParkerLab/Mitchell_2013_ParkerLab.xlsx")
RBP1010=RBP10[,c(2,1)]
colnames(RBP1010)=c("Gene","ENSEMBL")
RBP1010$Gene=toupper(RBP1010$Gene)

#Generate comparison matrix for these
RBPlistGene = list(as.character(RBP11$Gene), as.character(RBP22$Gene), as.character(RBP33$Gene), as.character(RBP44$Gene), as.character(RBPKlass$Gene), as.character(RBP88$Gene), as.character(RBP99$Gene), as.character(RBP1010$Gene)) #List of all the gene names
RBPlistENSEMBL = list(as.character(RBP11$ENSEMBL), as.character(RBP22$ENSEMBL), as.character(RBP33$ENSEMBL), as.character(RBP44$ENSEMBL), as.character(RBPKlass$ENSEMBL), as.character(RBP88$ENSEMBL), as.character(RBP99$ENSEMBL), as.character(RBP1010$ENSEMBL))
RBPlistGene=setNames(RBPlistGene,c("Beckman 2017","Brannan 2016","Gonzalez 2015","Hogan 2008","Klass 2013","Scherrer 2010","Tsvetanova 2010","Mitchell 2013"))
RBPlistENSEMBL=setNames(RBPlistENSEMBL,c("Beckman 2017","Brannan 2016","Gonzalez 2015","Hogan 2008","Klass 2013","Scherrer 2010","Tsvetanova 2010","Mitchell 2013"))
matrixGene = make_comb_mat(RBPlistGene,mode = "distinct")
matrixENSEMBL = make_comb_mat(RBPlistENSEMBL) #These can be used to upset plot if needed
#Do it on both Gene and ENSEMBLE as sometimes the data is missing in one of the columns. Then at the end we can Rbind based on Gene or ENSEMBLE (make sure to keep all) to get all of the proteins back
simpleMatGene = list_to_matrix(RBPlistGene) %>%
  as.data.frame() #Comparison matrix that is readable and subsettable
simpleMatENSEMBLE = list_to_matrix(RBPlistENSEMBL) %>%
  as.data.frame()
simpleMatGene$Sum = rowSums(simpleMatGene)
simpleMatENSEMBLE$Sum = rowSums(simpleMatENSEMBLE) #Sum the rows - if the sum is 2 or more, it has appeared in atleast 2 datasets. We can the filter based on confidence, dependent upon the number of datasets that particular protein has been observed in
conf2RBPGene = which(simpleMatGene$Sum>1) %>%
  rownames(simpleMatGene)[.]
conf2RBPENSEMBLE = which(simpleMatENSEMBLE$Sum>1) %>%
  rownames(simpleMatENSEMBLE)[.]
ENSEMBLE=c()
Gene=c()
for(i in 1:length(conf2RBPENSEMBLE)){
  ENSEMBLE[i] = conf2RBPENSEMBLE[i]
  Gene[i] = as.data.frame(org.Sc.sgdCOMMON2ORF)[min(which(as.data.frame(org.Sc.sgdCOMMON2ORF)$systematic_name == conf2RBPENSEMBLE[i])),1]
  ENSEMBLEDF=data.frame(ENSEMBLE=ENSEMBLE,Gene=Gene)
}
ENSEMBLE=c()
Gene=c()
for(i in 1:length(conf2RBPGene)){
  Gene[i] = conf2RBPGene[i]
  ENSEMBLE[i] = as.data.frame(org.Sc.sgdCOMMON2ORF)[min(which(as.data.frame(org.Sc.sgdCOMMON2ORF)$gene_name == conf2RBPGene[i])),2]
  GENEDF=data.frame(ENSEMBLE=ENSEMBLE,Gene=Gene)
}
conf2RBPS = merge(GENEDF,ENSEMBLEDF,by="ENSEMBLE",all=TRUE)
write.csv(conf2RBPS,"~/SGA_BioInfo/conf2RBPs.csv")

##############Load in SGA screen data ############
SGA2520 = read.xls("/Users/mfbx2cb8/Google Drive/SGA screen/CB/2520_all.xlsx")
sga2520trim = SGA2520[,c(1,5)]
for(i in 1:length(sga2520trim$ORF.Name)){
  sga2520trim$GENE[i] = as.data.frame(org.Sc.sgdCOMMON2ORF)[min(which(as.data.frame(org.Sc.sgdCOMMON2ORF)$systematic_name == as.character(sga2520trim$ORF.Name[i]))),1]
}
colnames(sga2520trim) = c("ORF.Name", "Gene","GENE")
SGA2521 = read.xls("/Users/mfbx2cb8/Google Drive/SGA screen/CB/2521_all.xlsx")
sga2521trim = SGA2521[,c(1,5)]
for(i in 1:length(sga2521trim$ORF.name)){
  sga2521trim$GENE[i] = as.data.frame(org.Sc.sgdCOMMON2ORF)[min(which(as.data.frame(org.Sc.sgdCOMMON2ORF)$systematic_name == as.character(sga2521trim$ORF.name[i]))),1]
}
colnames(sga2521trim) = c("ORF.Name", "Gene","GENE")
#It seems that some of the IDs and common names don't line up correctly in the datasets, so all analysis will be perofrmed on common name and ID
#Generate a list of those that impact both granule types, this will have repeated instances in it, but we can sort that out later.

sga25212520 = merge(sga2521trim,sga2520trim,by="ORF.Name")

######Compare SGA to RBPs #######
sga_both_RBPENSEMBLE=sga25212520[which(as.character(sga25212520$ORF.Name) %in% conf2RBPS$ENSEMBLE),]#Identified in both screens and are RBPs
sga_2520_RBP=sga2520trim[which(as.character(sga2520trim$ORF.Name) %in% conf2RBPS$ENSEMBLE),]
sga_2521_RBP=sga2521trim[which(as.character(sga2521trim$ORF.Name) %in% conf2RBPS$ENSEMBLE),]
sga_2520_RBP$SGD=bitr(sga_2520_RBP$GENE, fromType = "GENENAME", toType = c("SGD"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,2]
sga_2521_RBP$SGD=bitr(sga_2521_RBP$GENE, fromType = "GENENAME", toType = c("SGD"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,2]
sga_both_RBPENSEMBLE$SGD=bitr(sga_both_RBPENSEMBLE$GENE.x, fromType = "GENENAME", toType = c("SGD"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,2]
sga_2520_RBP$ENTREZ=bitr(sga_2520_RBP$GENE, fromType = "GENENAME", toType = c("ENTREZID"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,2]
sga_2521_RBP$ENTREZ=bitr(sga_2521_RBP$GENE, fromType = "GENENAME", toType = c("ENTREZID"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,2]
sga_both_RBPENSEMBLE$ENTREZ=bitr(sga_both_RBPENSEMBLE$GENE.x, fromType = "GENENAME", toType = c("ENTREZID"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,2]
sga_2520_RBP$ENSEMBL=bitr(sga_2520_RBP$GENE, fromType = "GENENAME", toType = c("ENSEMBL"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,2]
sga_2521_RBP$ENSEMBL=bitr(sga_2521_RBP$GENE, fromType = "GENENAME", toType = c("ENTREZID"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,2]
sga_both_RBPENSEMBLE$ENSEMBL=bitr(sga_both_RBPENSEMBLE$GENE.x, fromType = "GENENAME", toType = c("ENTREZID"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,2]
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/NIP1RBPs.csv", sga_2520_RBP)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/PDC1RBPs.csv", sga_2521_RBP)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/CommonRBPs.csv", sga_both_RBPENSEMBLE)

ego2520RBP_BP = enrichGO(gene=sga_2520_RBP$GENE, OrgDb = org.Sc.sgd.db, ont="BP", keyType = "GENENAME", pAdjustMethod = "bonferroni",pool = TRUE)
ego2520RBP_MF =  enrichGO(gene=sga_2520_RBP$GENE, OrgDb = org.Sc.sgd.db, ont="MF", keyType = "GENENAME",pAdjustMethod = "bonferroni",pool = TRUE)
ego2520RBP_CC =  enrichGO(gene=sga_2520_RBP$GENE, OrgDb = org.Sc.sgd.db, ont="CC", keyType = "GENENAME",pAdjustMethod = "bonferroni",pool = TRUE)
ego2521RBP_BP = enrichGO(gene=sga_2521_RBP$GENE, OrgDb = org.Sc.sgd.db, ont="BP", keyType = "GENENAME",pAdjustMethod = "bonferroni",pool = TRUE)
ego2521RBP_MF =  enrichGO(gene=sga_2521_RBP$GENE, OrgDb = org.Sc.sgd.db, ont="MF", keyType = "GENENAME",pAdjustMethod = "bonferroni",pool = TRUE)
ego2521RBP_CC =  enrichGO(gene=sga_2521_RBP$GENE, OrgDb = org.Sc.sgd.db, ont="CC", keyType = "GENENAME",pAdjustMethod = "bonferroni",pool = TRUE)
egobothRBP_BP = enrichGO(gene=sga_both_RBPENSEMBLE$GENE.x, OrgDb = org.Sc.sgd.db, ont="BP", keyType = "GENENAME",pAdjustMethod = "bonferroni",pool = TRUE)
egobothRBP_MF =  enrichGO(gene=sga_both_RBPENSEMBLE$GENE.x, OrgDb = org.Sc.sgd.db, ont="MF", keyType = "GENENAME",pAdjustMethod = "bonferroni",pool = TRUE)
egobothRBP_CC =  enrichGO(gene=sga_both_RBPENSEMBLE$GENE.x, OrgDb = org.Sc.sgd.db, ont="CC", keyType = "GENENAME",pAdjustMethod = "bonferroni",pool = TRUE)

write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/NIP1RBPGO_BP.csv", ego2520RBP_BP)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/NIP1RBPGO_MF.csv", ego2520RBP_MF)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/NIP1RBPGO_CC.csv", ego2520RBP_CC)

write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/PDC1RBPGO_BP.csv", ego2521RBP_BP)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/PDC1RBPGO_MF.csv", ego2521RBP_MF)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/PDC1RBPGO_CC.csv", ego2521RBP_CC)

write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/IntersectRBPGO_BP.csv", egobothRBP_BP)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/IntersectRBPGO_MF.csv", egobothRBP_MF)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/IntersectRBPGO_CC.csv", egobothRBP_CC)

ego_RBP_Combined_CC= merge_result(list(NIP=ego2520RBP_CC, PDC=ego2521RBP_CC, intersect=egobothRBP_CC))

p_ego_RBP_Combined_CC = ggplot(subset(ego_RBP_Combined_CC[], pvalue<0.01) , aes(x=Count,y=reorder(Description,Count)))+
  geom_point(shape=21, aes(fill=-log10(qvalue),size=-log10(p.adjust)))+
  theme_Publication_GridlinesTesting()+
  facet_wrap(~Cluster)+
  labs(size=NULL)+
  theme(axis.text.y = element_text(size=5))+
  scale_fill_gradient(low = "blue", high="red")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits=c(0,10),breaks=c(seq(0,10,2)))

dev.off()
cairo_pdf("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/GO_CC.pdf", width=5, height=4)
p_ego_RBP_Combined_CC
dev.off()

vego_RBP_Combined_BP= merge_result(list(NIP=ego2520RBP_BP, PDC=ego2521RBP_BP, intersect=egobothRBP_BP))

p_ego_RBP_Combined_BP = ggplot(subset(ego_RBP_Combined_BP[], pvalue<0.01) , aes(x=Count,y=reorder(Description,Count)))+
  geom_point(shape=21, aes(fill=-log10(qvalue),size=p.adjust))+
  theme_Publication_GridlinesTesting()+
  facet_wrap(~Cluster)+
  labs(size=NULL)+
  theme(axis.text.y = element_text(size=5))+
  scale_fill_gradient(low = "blue", high="red")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))


dev.off()
cairo_pdf("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/GO_CC.pdf", width=7.5, height=7.5)
p_ego_Combined
dev.off()
#pHyper analysis for overenrichment of RBPs
m=length(conf2RBPS[,1]) #Number of marks, e.g. the total number of genes with RBP function annotated
N=6275 #Total number of genes in yeast (in an ideal world, this number would be less, e.g. the number of proteins potentially discoverable by the methods)
n=N-m #non-marked elements
k=length(sga2520trim[,1]) #Background number of genes we were analysing - e.g. the number of genes identified in 2520
x=length(sga_2520_RBP[,1]) #Number of marked genes in background
p.value2520=phyper(q=x-1,m=m,n=n,k=k,lower.tail = FALSE) #Do x-1 because want to calculate p-val for enrichment. Therefore if p is signif for x-1, we are enriched.

k=length(sga2521trim[,1]) #Background number of genes we were analysing - e.g. the number of genes identified in 2520
x=length(sga_2521_RBP[,1]) #Number of marked genes in background
p.value2521=phyper(q=x-1,m=m,n=n,k=k,lower.tail = FALSE) #Do x-1 because want to calculate p-val for enrichment. Therefore if p is signif for x-1, we are enriched.

k=length(sga25212520[,1]) #Background number of genes we were analysing - e.g. the number of genes identified in 2520
x=length(sga_both_RBPENSEMBLE[,1]) #Number of marked genes in background
p.valueBoth=phyper(q=x-1,m=m,n=n,k=k,lower.tail = FALSE) #Do x-1 because want to calculate p-val for enrichment. Therefore if p is signif for x-1, we are enriched.

RBPEnrichment = data.frame(sample=c("2520","2521","both","background"), 
                           Total=c(length(sga2520trim[,1]),length(sga2521trim[,1]),length(sga25212520[,1]),N),
                           RBPs=c(length(sga_2520_RBP[,1]),length(sga_2521_RBP[,1]),length(sga_both_RBPENSEMBLE[,1]),m),
                           prop=c(length(sga_2520_RBP[,1])/length(sga2520trim[,1]),length(sga_2521_RBP[,1])/length(sga2521trim[,1]),length(sga_both_RBPENSEMBLE[,1])/length(sga25212520[,1]),m/N),
                           pval=c(p.value2520,p.value2521,p.valueBoth,"NA"))

RBPEnrichment$sample = factor(as.factor(RBPEnrichment$sample), levels=c("background","2520","2521","both"))
ggplot(RBPEnrichment,aes(x=sample, y=prop))+
  geom_bar(stat="identity")+
  theme_Publication_Gridlines()+
  scale_y_continuous(limits = c(0,1))

#####RBP PANTHER STUFF#####
NIP1_Panther = read.table("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/Panther_NIP1/2520pantherGeneList.txt",header=F,fill = TRUE, sep="\t")
PDC1_Panther = read.table("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/Panther_PDC1/pantherGeneList2521.txt",header=F,fill = TRUE, sep="\t")
Common_Panther = read.table("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/Panther_Common/pantherGeneList.txt",header=F,fill = TRUE, sep="\t")
NIP1_Panther$Cond = "NIP1"
PDC1_Panther$Cond =  "PDC1"
#Common missing V7 for some reason - just add it its filled with nas anyway
Common_Panther$V7 = NA
Common_Panther$Cond = "Common"
PantherLists_Collated = rbind(NIP1_Panther, PDC1_Panther, Common_Panther)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/PantherCollatedList.csv", PantherLists_Collated)

NIP1Chart = read.table("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/Panther_NIP1/2520pantherChart.txt",header=F,fill = TRUE, sep="\t")
PDCChart = read.table("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/Panther_PDC1/pantherChart.txt",header=F,fill = TRUE, sep="\t")
CommonChart = read.table("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/Panther_Common/pantherChart.txt",header=F,fill = TRUE, sep="\t")

NIP1Chart$Cond = "NIP1" 
PDCChart$Cond = "PDC1"
CommonChart$Cond = "Intersect" 

pantherCollated = rbind(NIP1Chart, PDCChart, CommonChart)
pantherCollated$V5 = pantherCollated$V5 %>% substr(., start = 1,stop=nchar(as.character(.))-1) %>% as.numeric()
pantherCollated$V4 = pantherCollated$V4 %>% substr(., start = 1,stop=nchar(as.character(.))-1) %>% as.numeric()
pantherCollated$Cond = factor(as.factor(pantherCollated$Cond),levels=c("NIP1","PDC1","Intersect"))

write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/PantherCollatedList.csv", PantherLists_Collated)

PantherPlot = ggplot(pantherCollated, aes(x = Cond, y = V5, fill = V2)) + 
  geom_bar(position = "fill",stat = "identity", col="black") +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#FFFFFF"))+
  theme_Publication_Gridlines(base_family = "Avenir Black",legendPos = "right",legendDir = "vertical")+
  labs(x="Dataset",y="Percentage of dataset",fill = "Panther class")

dev.off()
cairo_pdf("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/PantherPlot.pdf", width = 6.33, height = 4.77)
PantherPlot
dev.off()


#####RBP Disorder#####
NIP1_RBP_IDR = subset(sga_2520_RBP, sga_2520_RBP$ORF.Name %in% IDRGenesStrict$Seqid)
PDC1_RBP_IDR = subset(sga_2521_RBP, sga_2521_RBP$ORF.Name %in% IDRGenesStrict$Seqid)

NIP1_NotTranslation_NotNucleicAcidMetabolism = subset(sga_2520_RBP, ORF.Name %in% c("YGL008C","YER133W","YPR183W","YDR224C","YFR005C","YNL064C","YMR186W"))
PDC1_NotTranslation_NotNucleicAcidMetabolism = subset(sga_2521_RBP, ORF.Name %in% c("YER155C","YBR196C","YHR183W","YOR216C","YNL064C"))

NIP1_NotTranslation_NotNucleicAcidMetabolism_IDR = subset(NIP1_NotTranslation_NotNucleicAcidMetabolism, ORF.Name %in% IDRGenesStrict$Seqid)
PDC1_NotTranslation_NotNucleicAcidMetabolism_IDR = subset(PDC1_NotTranslation_NotNucleicAcidMetabolism, ORF.Name %in% IDRGenesStrict$Seqid)

write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/NIP1_RBPs_IDRs.csv", NIP1_RBP_IDR)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/RBPs/PDC1_RBPs_IDRs.csv", PDC1_RBP_IDR)

draw.pairwise.venn(area1 = length(sga_2520_RBP$Gene),
                   area2 = length(IDRGenesStrict$Seqid),
                   cross.area = length(intersect(sga_2520_RBP$ORF.Name,IDRGenesStrict$Seqid)))

draw.pairwise.venn(area1 = length(sga_2521_RBP$Gene),
                   area2 = length(IDRGenesStrict$Seqid),
                   cross.area = length(intersect(sga_2521_RBP$ORF.Name,IDRGenesStrict$Seqid)))

draw.triple.venn(area1 = length(sga_2520_RBP$Gene),
                 area2 = length(sga_2521_RBP$Gene),
                 area3 = length(IDRGenesStrict$Seqid),
                 n12 = length(intersect(sga_2520_RBP$ORF.Name,sga_2521_RBP$ORF.Name)),
                 n13 = length(intersect(sga_2520_RBP$ORF.Name,IDRGenesStrict$Seqid)),
                 n23 = length(intersect(sga_2521_RBP$ORF.Name,IDRGenesStrict$Seqid)),
                 n123 = length(Reduce(intersect, list(sga_2521_RBP$ORF.Name,sga_2520_RBP$ORF.Name,IDRGenesStrict$Seqid))))

########### Data from SGD ###########
#This data has information on the proteins encoded by these genes, such as their Gravy score, aliphatic score, abundance, length etc
#Read the data in
SGD2520 = read.csv("~/SGA_BioInfo/2520_LoadsOfInfo_Headed.csv", header=T)
SGD2521 = read.csv("~/SGA_BioInfo/2521_LoadsOfInfo_Headed.csv", header=T)
SGD = read.csv("~/SGA_BioInfo/YeastGenome_LoadsOfInfo_Headed.csv", header=T)  #Information on the whole yeast proteome - including dubious ORFs etc.
SGD_Common = SGD2521[which(SGD2521$Gene...Systematic.Name %in% SGD2520$Gene...Systematic.Name), ] #This is the intersect as determined by this function and the SGD website, although it only returns 12 overlaps, whereas the intersect we did above (sga25202521) has 45 overlaps.
SGD_Intersect = read.csv("~/SGA_BioInfo/Intersect_LoadsOfInfo_Headed.csv", header=T)
#Add an extra column that is a kind of 'id' based on where the data comes from - e.g. 2520, 2510, genome etc. Then we can just rbind them all togehter
SGD2520$ID = "NIP1"
SGD2521$ID = "PDC1"
SGD$ID = "BKGD"
SGD_Common$ID = "Int"
SGD_Intersect$ID = "Int"
SGD_All = rbind(SGD,SGD2520, SGD2521, SGD_Common)
SGD_All$ID = factor(as.factor(SGD_All$ID), levels=c("BKGD", "NIP1","PDC1","Int"))
p1=ggplot(SGD_All, aes(x=SGD_All$ID, y = SGD_All$Gene...Proteins...Codon.Bias, label=SGD_All$Gene...Standard.Name)) + geom_quasirandom(aes(col=ID))+geom_boxplot(fill=NA,col="black", notchwidth = 0.2, outlier.alpha=0) + theme_Publication_GridlinesTesting() + scale_color_manual(values=c("#845EC2","#D65DB1","#FF6F91","#939598"))+labs(y="Codon Bias", x=NULL)+theme(legend.position="none")
ggsave("~/SGA_BioInfo/Figs/Codon_bias.pdf",device = cairo_pdf)
p2=ggplot(SGD_All, aes(x=SGD_All$ID, y = SGD_All$Gene...Proteins...Molecular.Weight, label=SGD_All$Gene...Standard.Name)) + geom_quasirandom(aes(col=ID))+geom_boxplot(fill=NA,col="black", notchwidth = 0.2, outlier.alpha=0) + theme_Publication_GridlinesTesting() + scale_color_manual(values=c("#845EC2","#D65DB1","#FF6F91","#939598"))+labs(y="Molecular Weight (Da)", x=NULL)+theme(legend.position="none")
ggsave("~/SGA_BioInfo/Figs/MW.pdf",device = cairo_pdf)
p3=ggplot(SGD_All, aes(x=SGD_All$ID, y = SGD_All$Gene...Proteins...Protein.Abundance.Median.Absolute.Deviation, label=SGD_All$Gene...Standard.Name)) + geom_quasirandom(aes(col=ID))+geom_boxplot(fill=NA,col="black", notchwidth = 0.2, outlier.alpha=0) + theme_Publication_GridlinesTesting() + scale_color_manual(values=c("#845EC2","#D65DB1","#FF6F91","#939598"))+labs(y="Abundance (CPC)", x=NULL)+theme(legend.position="none")
ggsave("~/SGA_BioInfo/Figs/Abundance.pdf",device = cairo_pdf)
p4=ggplot(SGD_All, aes(x=SGD_All$ID, y = SGD_All$Gene...Proteins...Protein.Abundance.Median.Absolute.Deviation, label=SGD_All$Gene...Standard.Name)) + geom_quasirandom(aes(col=ID))+geom_boxplot(fill=NA,col="black", notchwidth = 0.2, outlier.alpha=0) + theme_Publication_GridlinesTesting() + scale_color_manual(values=c("#845EC2","#D65DB1","#FF6F91","#939598"))+scale_y_continuous(limits = c(0,0.1e+05))+labs(y="Abundance (CPC)", x=NULL)+theme(legend.position="none")
ggsave("~/SGA_BioInfo/Figs/AbundanceZoom.pdf",device = cairo_pdf)
p5=ggplot(SGD_All, aes(x=SGD_All$ID, y = SGD_All$Gene...Proteins...Gravy.Score, label=SGD_All$Gene...Standard.Name)) + geom_quasirandom(aes(col=ID))+geom_boxplot(fill=NA,col="black", notchwidth = 0.2, outlier.alpha=0) + theme_Publication_GridlinesTesting() + scale_color_manual(values=c("#845EC2","#D65DB1","#FF6F91","#939598"))+labs(y="Gravy Score", x=NULL)+theme(legend.position="none")
ggsave("~/SGA_BioInfo/Figs/Gravy.pdf",device = cairo_pdf)
p6=ggplot(SGD_All, aes(x=SGD_All$ID, y = SGD_All$Gene...Proteins...P.I, label=SGD_All$Gene...Standard.Name)) + geom_quasirandom(aes(col=ID))+geom_boxplot(fill=NA,col="black", notchwidth = 0.2, outlier.alpha=0) + theme_Publication_GridlinesTesting() + scale_color_manual(values=c("#845EC2","#D65DB1","#FF6F91","#939598"))+labs(y="Isoelectric Point", x=NULL)+theme(legend.position="none")
p7=ggplot(SGD_All, aes(x=SGD_All$ID, y = SGD_All$Gene...Proteins...Instability.Index, label=SGD_All$Gene...Standard.Name)) + geom_quasirandom(aes(col=ID))+geom_boxplot(fill=NA,col="black", notchwidth = 0.2, outlier.alpha=0) + theme_Publication_GridlinesTesting() + scale_color_manual(values=c("#845EC2","#D65DB1","#FF6F91","#939598"))+labs(y="Instability Index", x=NULL)+theme(legend.position="none")
ggsave("~/SGA_BioInfo/Figs/Instability.pdf",device = cairo_pdf)
p8=ggplot(SGD_All, aes(x=SGD_All$ID, y = SGD_All$Gene...Proteins...Aliphatic.Index, label=SGD_All$Gene...Standard.Name)) + geom_quasirandom(aes(col=ID))+geom_boxplot(fill=NA,col="black", notchwidth = 0.2, outlier.alpha=0) + theme_Publication_GridlinesTesting() + scale_color_manual(values=c("#845EC2","#D65DB1","#FF6F91","#939598"))+labs(y="Aliphatic Index", x=NULL)+theme(legend.position="none")
ggsave("~/SGA_BioInfo/Figs/Aliphatic.pdf",device = cairo_pdf)
p1+p2+p4+p5+p6+p7+plot_layout(ncol=3)& theme(plot.margin =margin(2,2,2,10))
ggsave("~/SGA_BioInfo/Figs/Fig.pdf",device = cairo_pdf,width = 160,height = 175,units = "mm")

#Translation related stuff
weinbergTranslation = xlsx::read.xlsx("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Translation/Weinberg Cell Reports 2016 Flash freeze.xlsx", sheetIndex =1)
weinbergTranslation$CommonName = bitr(weinbergTranslation$ORF, fromType = "ORF", toType = c("GENENAME"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,3]
weinbergTranslation2520 = subset(weinbergTranslation, weinbergTranslation$ORF %in% SGD2520$Gene...Systematic.Name)
weinbergTranslation2521 = subset(weinbergTranslation, weinbergTranslation$ORF %in% SGD2521$Gene...Systematic.Name)
weinbergTranslationCommon = subset(weinbergTranslation2521, weinbergTranslation2521$ORF %in% weinbergTranslation2520$ORF)
weinbergTranslation$ID = "BKGD"
weinbergTranslation2520$ID = "NIP1" 
weinbergTranslation2521$ID = "PDC1"
weinbergTranslationCommon$ID = "Int"
weinbergTranslation_All = rbind(weinbergTranslation,weinbergTranslation2520, weinbergTranslation2521, weinbergTranslationCommon)
weinbergTranslation_All$ID = factor(as.factor(weinbergTranslation_All$ID), levels=c("BKGD", "NIP1","PDC1","Int"))

#Plots for thesis
pMolWeight = ggplot(SGD_All, aes(x=ID, y = Gene...Proteins...Molecular.Weight, label=Gene...Standard.Name))+
  geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_GridlinesTesting(base_family = "Avenir Black")+
  scale_fill_manual(values=c("#C2C3C6","#F16058","#F4B368","#3F59A8")) +
  scale_colour_manual(values=c("#939598","#B84845","#D88D42","#263976")) +
  labs(y="Molecular weight (Da)", x="SGA hits")+
  theme(legend.position="none")+
  theme(plot.margin = unit(c(3, 3, 3, 3), "mm"))

pGravy = ggplot(SGD_All, aes(x=ID, y = Gene...Proteins...Gravy.Score, label=Gene...Standard.Name)) + 
  geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_GridlinesTesting(base_family = "Avenir Black")+
  scale_fill_manual(values=c("#C2C3C6","#F16058","#F4B368","#3F59A8")) +
  scale_colour_manual(values=c("#939598","#B84845","#D88D42","#263976")) +
  labs(y="Hydrophobicity (GRAVY Score)", x="SGA hits")+theme(legend.position="none")+
  theme(plot.margin = unit(c(3, 3, 3, 3), "mm"))


pIsoElec = ggplot(SGD_All, aes(x=ID, y = Gene...Proteins...P.I, label=Gene...Standard.Name)) + 
  geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_GridlinesTesting(base_family = "Avenir Black")+
  scale_fill_manual(values=c("#C2C3C6","#F16058","#F4B368","#3F59A8")) +
  scale_colour_manual(values=c("#939598","#B84845","#D88D42","#263976")) +
  labs(y="Isoelectric point (pI)", x="SGA hits")+theme(legend.position="none")+
  theme(plot.margin = unit(c(3, 3, 3, 3), "mm"))

pAbundance = ggplot(SGD_All, aes(x=ID, y = log10(Gene...Proteins...Protein.Abundance.Median.Absolute.Deviation), label=Gene...Standard.Name)) + 
  geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_GridlinesTesting(base_family = "Avenir Black")+
  scale_fill_manual(values=c("#C2C3C6","#F16058","#F4B368","#3F59A8")) +
  scale_colour_manual(values=c("#939598","#B84845","#D88D42","#263976")) +
  labs(y=expression("Protein abundance (Log"[10]*" molecules per cell)"), x="SGA hits")+theme(legend.position="none")+
  theme(plot.margin = unit(c(3, 3, 3, 3), "mm"))

pTE = ggplot(weinbergTranslation_All, aes(x=ID, y = TE..RPF.mRNA., label=ORF)) + 
  geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_GridlinesTesting(base_family = "Avenir Black")+
  scale_fill_manual(values=c("#C2C3C6","#F16058","#F4B368","#3F59A8")) +
  scale_colour_manual(values=c("#939598","#B84845","#D88D42","#263976")) +
  labs(y=expression("Ribosome occupancy (RPKM"[RPF]*"/RPKM"[mRNA]*")"), x="SGA hits")+theme(legend.position="none")+
  theme(plot.margin = unit(c(3, 3, 3, 3), "mm"))

pSynthesisRate = ggplot(weinbergTranslation_All, aes(x=ID, y = log10(SynthesisRate), label=ORF)) + 
  geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_GridlinesTesting(base_family = "Avenir Black")+
  scale_fill_manual(values=c("#C2C3C6","#F16058","#F4B368","#3F59A8")) +
  scale_colour_manual(values=c("#939598","#B84845","#D88D42","#263976")) +
  labs(y=expression("Log"[10]*" Protein synthesis rate "), x="SGA hits")+theme(legend.position="none")+
  theme(plot.margin = unit(c(3, 3, 3, 3), "mm"))

p1 <- ggplot_build(pMolWeight)
#Now edit the plot data - the [[x]] accession is the order of plotting, so if you plot the beeswarm 1st or 3rd it would be [[1]] or [[3]], respectively
#We need to remove points that fall to the left of the centre line. Again, this depends on the number of x-axis factors you have.
#The below function takes points which fall to the left of the x-label (e.g. x-group is -ve) and makes these positive.
#It leaves values which have x-group as +ve because these don't need to be changed
p1$data[[1]]=p1$data[[1]] %>% mutate(x=case_when(
  x-round(x) < 0 ~ abs(x-round(x))+round(x),
  TRUE ~ x)
)

p2 <- ggplot_build(pGravy)
#Now edit the plot data - the [[x]] accession is the order of plotting, so if you plot the beeswarm 1st or 3rd it would be [[1]] or [[3]], respectively
#We need to remove points that fall to the left of the centre line. Again, this depends on the number of x-axis factors you have.
#The below function takes points which fall to the left of the x-label (e.g. x-group is -ve) and makes these positive.
#It leaves values which have x-group as +ve because these don't need to be changed
p2$data[[1]]=p2$data[[1]] %>% mutate(x=case_when(
  x-round(x) < 0 ~ abs(x-round(x))+round(x),
  TRUE ~ x)
)

p3 <- ggplot_build(pIsoElec)
#Now edit the plot data - the [[x]] accession is the order of plotting, so if you plot the beeswarm 1st or 3rd it would be [[1]] or [[3]], respectively
#We need to remove points that fall to the left of the centre line. Again, this depends on the number of x-axis factors you have.
#The below function takes points which fall to the left of the x-label (e.g. x-group is -ve) and makes these positive.
#It leaves values which have x-group as +ve because these don't need to be changed
p3$data[[1]]=p3$data[[1]] %>% mutate(x=case_when(
  x-round(x) < 0 ~ abs(x-round(x))+round(x),
  TRUE ~ x)
)

p4 <- ggplot_build(pAbundance)
#Now edit the plot data - the [[x]] accession is the order of plotting, so if you plot the beeswarm 1st or 3rd it would be [[1]] or [[3]], respectively
#We need to remove points that fall to the left of the centre line. Again, this depends on the number of x-axis factors you have.
#The below function takes points which fall to the left of the x-label (e.g. x-group is -ve) and makes these positive.
#It leaves values which have x-group as +ve because these don't need to be changed
p4$data[[1]]=p4$data[[1]] %>% mutate(x=case_when(
  x-round(x) < 0 ~ abs(x-round(x))+round(x),
  TRUE ~ x)
)

p5 <- ggplot_build(pTE)
#Now edit the plot data - the [[x]] accession is the order of plotting, so if you plot the beeswarm 1st or 3rd it would be [[1]] or [[3]], respectively
#We need to remove points that fall to the left of the centre line. Again, this depends on the number of x-axis factors you have.
#The below function takes points which fall to the left of the x-label (e.g. x-group is -ve) and makes these positive.
#It leaves values which have x-group as +ve because these don't need to be changed
p5$data[[1]]=p5$data[[1]] %>% mutate(x=case_when(
  x-round(x) < 0 ~ abs(x-round(x))+round(x),
  TRUE ~ x)
)

p6 <- ggplot_build(pSynthesisRate)
#Now edit the plot data - the [[x]] accession is the order of plotting, so if you plot the beeswarm 1st or 3rd it would be [[1]] or [[3]], respectively
#We need to remove points that fall to the left of the centre line. Again, this depends on the number of x-axis factors you have.
#The below function takes points which fall to the left of the x-label (e.g. x-group is -ve) and makes these positive.
#It leaves values which have x-group as +ve because these don't need to be changed
p6$data[[1]]=p6$data[[1]] %>% mutate(x=case_when(
  x-round(x) < 0 ~ abs(x-round(x))+round(x),
  TRUE ~ x)
)


plot_grid(ggplot_gtable(p1),ggplot_gtable(p2),ggplot_gtable(p3),ggplot_gtable(p4),ggplot_gtable(p5),ggplot_gtable(p6),ncol=2,align="hv")



vt
####With labelled outliers#####
#Outlier function
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

weinbergTranslation_All%>%
  group_by(ID) %>%
  mutate(outlier=ifelse(is_outlier(log10(SynthesisRate)),CommonName, NA)) %>%
  ggplot( aes(x=ID, y = log10(SynthesisRate), label=ORF)) +
  geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  ggrepel::geom_text_repel(data=. %>% filter(!is.na(outlier) && ID != "BKGD"), aes(label=outlier))+
  geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_GridlinesTesting(base_family = "Avenir Black")+
  scale_fill_manual(values=c("#C2C3C6","#F16058","#F4B368","#3F59A8")) +
  scale_colour_manual(values=c("#939598","#B84845","#D88D42","#263976")) +
  labs(y=expression("Log"[10]*" Protein synthesis rate "), x="SGA hits")+theme(legend.position="none")+
  theme(plot.margin = unit(c(3, 3, 3, 3), "mm"))




#Save data into translation folder too!

write.csv(file = "~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Translation/2520_TranslationInfo.csv", weinbergTranslation2520)
write.csv(file = "~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Translation/2521_TranslationInfo.csv", weinbergTranslation2521)
write.csv(file = "~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Translation/Common_TranslationInfo.csv", weinbergTranslationCommon)

#Stat tests!

test1=wilcox.test(SynthesisRate~ID, data=subset(weinbergTranslation_All, ID=="BKGD"|ID=="NIP1"))
test2=wilcox.test(SynthesisRate~ID, data=subset(weinbergTranslation_All, ID=="BKGD"|ID=="PDC1"))
test3=wilcox.test(SynthesisRate~ID, data=subset(weinbergTranslation_All, ID=="BKGD"|ID=="Int"))
test4=wilcox.test(TE..RPF.mRNA.~ID, data=subset(weinbergTranslation_All, ID=="BKGD"|ID=="NIP1"))
test5=wilcox.test(TE..RPF.mRNA.~ID, data=subset(weinbergTranslation_All, ID=="BKGD"|ID=="PDC1"))
test6=wilcox.test(TE..RPF.mRNA.~ID, data=subset(weinbergTranslation_All, ID=="BKGD"|ID=="Int"))
test7=wilcox.test(Gene...Proteins...Molecular.Weight~ID, data=subset(SGD_All, ID=="BKGD"|ID=="NIP1"))
test8=wilcox.test(Gene...Proteins...Molecular.Weight~ID, data=subset(SGD_All, ID=="BKGD"|ID=="PDC1"))
test9=wilcox.test(Gene...Proteins...Molecular.Weight~ID, data=subset(SGD_All, ID=="BKGD"|ID=="Int"))
test10=wilcox.test(Gene...Proteins...Gravy.Score~ID, data=subset(SGD_All, ID=="BKGD"|ID=="NIP1"))
test11=wilcox.test(Gene...Proteins...Gravy.Score~ID, data=subset(SGD_All, ID=="BKGD"|ID=="PDC1"))
test12=wilcox.test(Gene...Proteins...Gravy.Score~ID, data=subset(SGD_All, ID=="BKGD"|ID=="Int"))
test13=wilcox.test(Gene...Proteins...P.I~ID, data=subset(SGD_All, ID=="BKGD"|ID=="NIP1"))
test14=wilcox.test(Gene...Proteins...P.I~ID, data=subset(SGD_All, ID=="BKGD"|ID=="PDC1"))
test15=wilcox.test(Gene...Proteins...P.I~ID, data=subset(SGD_All, ID=="BKGD"|ID=="Int"))
test16=wilcox.test(Gene...Proteins...Protein.Abundance.Median.Absolute.Deviation~ID, data=subset(SGD_All, ID=="BKGD"|ID=="NIP1"))
test17=wilcox.test(Gene...Proteins...Protein.Abundance.Median.Absolute.Deviation~ID, data=subset(SGD_All, ID=="BKGD"|ID=="PDC1"))
test18=wilcox.test(Gene...Proteins...Protein.Abundance.Median.Absolute.Deviation~ID, data=subset(SGD_All, ID=="BKGD"|ID=="Int"))

pValues = c(test1$p.value, test2$p.value, test3$p.value,test4$p.value, test5$p.value, test6$p.value,test7$p.value, test8$p.value, test9$p.value,test10$p.value, test11$p.value, test12$p.value,test13$p.value, test14$p.value, test15$p.value,test16$p.value, test17$p.value, test18$p.value)
test = c(rep("SynthesisRate",3),rep("TE",3),rep("MolWeight",3),rep("Gravy",3),rep("P.I",3),rep("Abundance",3))
ID = c(rep(c("NIP1","PDC1","Int"),6))
pValDF = data.frame(pValues=pValues, data=test, ID=ID)
#Have to arrange is as pAdj outputs in order  - below orders within each group
pValDF = pValDF %>% group_by(data) %>% arrange(pValues, .by_group = TRUE)

#Have to splot the p adjustments into each data measured - e.g. MW has its own p adjustments
pValDF1=p.adjust(pValDF$pValues[1:3],method = "bonferroni")
pValDF2=p.adjust(pValDF$pValues[4:6],method = "bonferroni")
pValDF3=p.adjust(pValDF$pValues[7:9],method = "bonferroni")
pValDF4=p.adjust(pValDF$pValues[10:12],method = "bonferroni")
pValDF5=p.adjust(pValDF$pValues[13:15],method = "bonferroni")
pValDF6=p.adjust(pValDF$pValues[16:18],method = "bonferroni")

pValDF$p.Adj=c(pValDF1,pValDF2,pValDF3,pValDF4,pValDF5,pValDF6)

write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/pValuesPhysicoChemicalData.csv", pValDF)



   
######GO Analysis using file made by Vicky - see her thesis for how she did it#######
GO_Screen=read.csv("~/SGA_BioInfo/GO_PizzingaAnalysis_Condensed.csv", header =T) #Read the file
GO_Screen = GO_Screen[with(GO_Screen, order(GO_Screen$p)), ] #order them by P-value
axisPalette = c("#C70039","#29335C")#Generate colour palette for x axis labels so we know which database they came from
names(axisPalette)=levels(GO_Screen$Database)#name the colours to attribute to the different databases
plot = ggplot(GO_Screen, aes(GO_Screen$mRNA,factor(GO_Screen$GO.Name, levels=unique(GO_Screen$GO.Name[order(GO_Screen$Database)]), ordered=TRUE)))
plot + geom_tile(aes(fill = -log10(GO_Screen$p)),colour = "black", size = 1, width =1, height = 1)+scale_fill_gradient(high='#C70039', low="#FFFFFF")+theme_Publication_GridlinesTesting()
plot + geom_point(aes(fill = -log10(GO_Screen$p),col=GO_Screen$Database), size = 2, shape=21)+scale_fill_gradient(high='#C70039', low="#FFFFFF")+theme_Publication_GridlinesTesting()+theme(axis.text = element_text(size = 4),axis.text.y = element_text(colour=axisPalette[GO_Screen$Database[order(GO_Screen$Database)]])) #This line sets the axis text colour
ggsave("~/SGA_BioInfo/Figs/GOAnalysis.pdf", device=cairo_pdf)

#####GO analysis myself using catSGA#####
#correctIDs<-bitr(Dcp1clusterMembers, fromType = "ORF", toType = c("GENENAME"), OrgDb = org.Sc.sgd.db)
#ego<-enrichGO(gene=correctIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont = "CC", keyType = "GENENAME", minGSSize=1)
#DOSE::enrichDO()
require(clusterProfiler)
require(org.Sc.sgd.db)
require(ggplot2)
require(scales)
library("rrvgo")
catSGA=readxl::read_xlsx("~/SGA_BioInfo/Categorised_Hits_Pos_Neg_Misloc_CB.xlsx")
SGAIDs = bitr(catSGA$SeqID, fromType = "ORF", toType = c("GENENAME"), OrgDb = org.Sc.sgd.db,)
egoSGA_BP = enrichGO(gene=SGAIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont="BP", keyType = "GENENAME",pAdjustMethod = "none",pool = TRUE)
egoSGA_MF = enrichGO(gene=SGAIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont="MF", keyType = "GENENAME",pAdjustMethod = "none",pool = TRUE)
egoSGA_CC = enrichGO(gene=SGAIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont="CC", keyType = "GENENAME",pAdjustMethod = "none",pool = TRUE)
NIPIDs = bitr(subset(catSGA, catSGA$mRNA=="NIP1")$SeqID, fromType = "ORF", toType = c("GENENAME"), OrgDb = org.Sc.sgd.db)
egoNIP_BP = enrichGO(gene=NIPIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont="BP", keyType = "GENENAME",pAdjustMethod = "none",pool = TRUE)
egoNIP_MF = enrichGO(gene=NIPIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont="MF", keyType = "GENENAME",pAdjustMethod = "none",pool = TRUE)
egoNIP_CC = enrichGO(gene=NIPIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont="CC", keyType = "GENENAME",pAdjustMethod = "none",pool = TRUE)
PDCIDs = bitr(subset(catSGA, catSGA$mRNA=="PDC1")$SeqID, fromType = "ORF", toType = c("GENENAME"), OrgDb = org.Sc.sgd.db)
egoPDC_BP = enrichGO(gene=PDCIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont="BP", keyType = "GENENAME",pAdjustMethod = "none",pool = TRUE)
egoPDC_MF = enrichGO(gene=PDCIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont="MF", keyType = "GENENAME",pAdjustMethod = "none",pool = TRUE)
egoPDC_CC = enrichGO(gene=PDCIDs$GENENAME, OrgDb = org.Sc.sgd.db, ont="CC", keyType = "GENENAME",pAdjustMethod = "none",pool = TRUE)

#Merge the MFs together for plotting
ego_Combined= merge_result(list(SGA=egoSGA_BP, NIP=egoNIP_BP, PDC=egoPDC_BP))

p_ego_Combined = ggplot(subset(ego_Combined[], pvalue<0.01) , aes(x=Count,y=reorder(Description,Count)))+
  geom_point(shape=21, aes(fill=-log10(qvalue)))+
  theme_Publication_GridlinesTesting()+
  facet_wrap(~Cluster)+
  labs(size=NULL)+
  theme(axis.text.y = element_text(size=5))+
  scale_fill_gradient(low = "blue", high="red")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))

dev.off()
cairo_pdf("~/Documents/Thesis--sssss/Chapter2/4.1 A high-throughput SGA screen/GO_Combined.pdf", width=7.5, height=7.5)
p_ego_Combined
dev.off()

semSim = GOSemSim::godata(org.Sc.sgd.db, ont="BP")

#Reduce number of ontologies by using semantic similarity with rrvgo
simMatrix_SGA_BP <- calculateSimMatrix(egoSGA_BP$ID,
                                       orgdb="org.Sc.sgd.db",
                                       ont="BP",
                                       method="Rel",
                                       semdata = semSim)
scores_SGA_BP <- setNames(-log10(egoSGA_BP$qvalue), egoSGA_BP$ID)
reducedTerms_SGA_BP <- reduceSimMatrix(simMatrix_SGA_BP,
                                scores_SGA_BP,
                                threshold=0.7,
                                orgdb="org.Sc.sgd.db")

simMatrix_NIP_BP <- calculateSimMatrix(egoNIP_BP$ID,
                                orgdb="org.Sc.sgd.db",
                                ont="BP",
                                method="Rel",
                                semdata = semSim)
scores_NIP_BP <- setNames(-log10(egoNIP_BP$qvalue), egoNIP_BP$ID)
reducedTerms_NIP_BP <- reduceSimMatrix(simMatrix_NIP_BP,
                                scores_NIP_BP,
                                threshold=0.7,
                                orgdb="org.Sc.sgd.db")

simMatrix_PDC_BP <- calculateSimMatrix(egoPDC_BP$ID,
                                       orgdb="org.Sc.sgd.db",
                                       ont="BP",
                                       method = "Wang",
                                       semdata = semSim)
scores_PDC_BP <- setNames(-log10(egoPDC_BP$qvalue), egoPDC_BP$ID)
reducedTerms_PDC_BP <- reduceSimMatrix(simMatrix_PDC_BP,
                                scores_PDC_BP,
                                threshold=0.7,
                                orgdb="org.Sc.sgd.db")
reducedTerms_SGA_BP$Cond = "Comb"
reducedTerms_NIP_BP$Cond = "NIP1"
reducedTerms_PDC_BP$Cond = "PDC1"

p1 = scatterPlot(simMatrix_SGA_BP,reducedTerms_SGA_BP,labelSize = 2, max.overlaps=) +theme_Publication_Gridlines()
p2 = scatterPlot(simMatrix_NIP_BP,reducedTerms_NIP_BP) +theme_Publication_Gridlines()
p3 = scatterPlot(simMatrix_PDC_BP,reducedTerms_PDC_BP) +theme_Publication_Gridlines()

reducedTerms_Combined = rbind(reducedTerms_SGA_BP,reducedTerms_NIP_BP, reducedTerms_PDC_BP)
reducedTerms_Combined$Cond = factor(as.factor(reducedTerms_Combined$Cond), levels=c("Comb","NIP1","PDC1"))
ggplot(reducedTerms_Combined, aes(x=Cond, y=parentTerm))+
  geom_point(shape=21, aes(size=score))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  theme(axis.text.y = element_text(size=5))+
  theme_Publication_Gridlines()

simplifiedReducedTerms_Combined = reducedTerms_Combined %>% group_by(Cond, parentTerm) %>% dplyr::summarise(score=sum(score),cluster=cluster)
#So will want plots pre-semantic reduction and following
cairo_pdf("~/Documents/Thesis--sssss/Chapter2/4.1 A high-throughput SGA screen/GO_SemSim_Summed.pdf", width=8, height =7.5)
ggplot(simplifiedReducedTerms_Combined, aes(x=Cond, y=reorder(parentTerm, cluster)))+
  geom_point(shape=21, aes(fill=Cond,size=score))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  theme(axis.text.y = element_text(size=5))+
  theme_Publication_Gridlines(base_family = "Avenir Black")+
  scale_fill_manual(values = c("#F1F0F0","#F16058","#F4B368"))+
  theme(axis.text.y = element_text(size=7),axis.text.x = element_text(size=7))

  
dev.off()


#########On combined dataset so all 3 are reduced in complexity together########
ego_Combined= merge_result(list(SGA=egoSGA_BP, NIP=egoNIP_BP, PDC=egoPDC_BP))
ego_Combined@compareClusterResult$Cluster=interaction(ego_Combined[]$Cluster,ego_Combined[]$qvalue)
simMatrix <- calculateSimMatrix(ego_Combined[]$ID,
                                orgdb="org.Sc.sgd.db",
                                ont="BP",
                                method="Rel",
                                semdata = semSim)
scores <- setNames(ego_Combined[]$Cluster, ego_Combined[]$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Sc.sgd.db")
reducedTerms=transform(reducedTerms, Cond = substr(score, 1, 3), score = as.numeric(substr(score, 5, nchar(as.character(score)))))
reducedTerms = reducedTerms %>% group_by(Cond, parentTerm) %>% dplyr::summarise(score=sum(-log10(score)),cluster=cluster, go=go)
reducedTerms$Cond = factor(as.factor(reducedTerms$Cond), levels=c("SGA","NIP","PDC"))
cairo_pdf("~/Documents/Thesis--sssss/Chapter2/4.1 A high-throughput SGA screen/GO_SemSim_AllTogether.pdf", width=5, height =5)
ggplot(reducedTerms, aes(x=Cond, y=reorder(parentTerm, go)))+
  geom_point(shape=21, aes(fill=Cond, size=as.numeric(score)))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  theme(axis.text.y = element_text(size=5))+
  theme_Publication_Gridlines(base_family = "Avenir Black")+
  scale_fill_manual(values = c("#F1F0F0","#F16058","#F4B368"))+
  theme_Publication_Gridlines()+
  theme(axis.text.y = element_text(size=7),axis.text.x = element_text(size=7))
dev.off()



#######Disorder Stuff#########
#Start with MobiDB - it was output as a .json, but i have converted it to CSV with an online tool. The output isn't
#ideal as it have loads of NULL values. Essentially, it lists the protein uniprot acc, then has a cell per amino acid, with a 
#Score for that amino acid in the 'mobidb_consensus_disorder_predictor_score' column. There is also a column with the disordered regions
#listed with the amino acid boundaries for these regions. e.g.  55 70 D means disordered region between aa 55 and 70. This
#is the 'mobidb_consensus_disorder_predictors_regions' column. The other columns seem pretty pointless.
#read in the csv (only 5000 rows to start off with so not working with a beast of a file)
disorder=read.csv("~/Downloads/Disorder Stuff/sqlify-result-b299b76c965a4.csv",header=T,stringsAsFactors = F)
disorderProcessed = list() #Store them in a list where each item is a protein - this will aid plotting indivdual profiles etc later :)
for(i in 1:length(disorder[,1])){
    if(i==1){
      x=0 #Set gene counter to 0
    }
    if(disorder[i,1]!="NULL"){  #If this is not NULL then it is a acc number and therefore the beginning of a new gene
      f=1 #f is going to be the counter for each gene
      k=1
      j=1
      x=x+1                     #So bump the counter by 1 to start a new item in the list
      print(x)
      disorderProcessed[x] = disorder[i,1]     #The name of this item is the acc, but not if its the first gene
      disorderProcessed[[x]]$Name=disorder[i,1]
      if(i<tail(which(disorder$acc!="NULL"),1)){
        disorderProcessed[[x]]$DisorderedRegions = rep(0,which(disorder$acc!="NULL")[x+1]-which(disorder$acc!="NULL")[x])  #Need fill in the blanks on the disorder thing.
      }
      if(i==tail(which(disorder$acc!="NULL"),1)){
        disorderProcessed[[x]]$DisorderedRegions = rep(0,length(disorder[,1])-tail(which(disorder$acc!="NULL"),1))  #Need to go back and modify the previous gene to fill in the blanks on the disorder thing.
      }
    }
  disorderProcessed[[x]]$AAScore[f] = as.numeric(disorder$mobidb_consensus_disorder_predictors_scores[i])     #The score is in the 6th column, could access it with a $ ident too.
  disorderProcessed[[x]]$AAScoreDerived[j] = as.numeric(disorder$mobidb_consensus_disorder_derived_scores_n[i])     #The score is in the 6th column, could access it with a $ ident too.
  if(!is.na(as.numeric(disorder$mobidb_consensus_disorder_predictors_regions[i])) & !is.na(as.numeric(disorder$mobidb_consensus_disorder_predictors_regions[i+1])) & disorder$mobidb_consensus_disorder_predictors_regions[i+2]=="D"){
    #If the value is not an NA, it must be a number, so we want to keep it. but only if the acc before it is not numeric! otherwise we'd end up with disordered regions trying to span regions between disordered regions (if that makes sense?)
    disorderProcessed[[x]]$DisorderedRegions[as.numeric(disorder$mobidb_consensus_disorder_predictors_regions[i]):as.numeric(disorder$mobidb_consensus_disorder_predictors_regions[i+1])] =1    #Set the values for these to 1 - a binary score of if it is disordered or not that can be overlayed with the amino acid seq & scores.
    disorderProcessed[[x]]$simpleRegions[k]=as.numeric(disorder$mobidb_consensus_disorder_predictors_regions[i])
    disorderProcessed[[x]]$simpleRegions[k+1]=as.numeric(disorder$mobidb_consensus_disorder_predictors_regions[i+1])
    k=k+2
  }
  if(!is.na(as.numeric(disorder$mobidb_consensus_disorder_derived_regions[i])) & !is.na(as.numeric(disorder$mobidb_consensus_disorder_derived_regions[i+1])) & disorder$mobidb_consensus_disorder_derived_regions[i+2]=="D"){
    #If the value is not an NA, it must be a number, so we want to keep it. but only if the acc before it is not numeric! otherwise we'd end up with disordered regions trying to span regions between disordered regions (if that makes sense?)
    disorderProcessed[[x]]$DisorderedRegionsStrict[as.numeric(disorder$mobidb_consensus_disorder_derived_regions[i]):as.numeric(disorder$mobidb_consensus_disorder_derived_regions[i+1])] =1    #Set the values for these to 1 - a binary score of if it is disordered or not that can be overlayed with the amino acid seq & scores.
    disorderProcessed[[x]]$simpleRegionsStrict[j]=as.numeric(disorder$mobidb_consensus_disorder_derived_regions[i])
    disorderProcessed[[x]]$simpleRegionsStrict[j+1]=as.numeric(disorder$mobidb_consensus_disorder_derived_regions[i+1])
    j=j+2
  }
  f=f+1
}

saveRDS(disorderProcessed,"~/SGA_BioInfo/DisorderProcessed")

SGD.sc=read.csv("~/SGA_BioInfo/org.sc.sgd_LoadsOfInfo_Headed.csv",header=T)
SGD.sc$ID="sc"
SGD_Combined = rbind(SGD,SGD.sc) 
SGD_Combined=SGD_Combined[!duplicated(SGD_Combined$Gene...Systematic.Name), ]
D2P2 = read.table("~/Downloads/Disorder Stuff/D2P2_Disorder_download.tsv",header=T,sep="\t",stringsAsFactors = F)
D2P2=D2P2[which(D2P2$Genome=="<i>Saccharomyces cerevisiae</i> <release>63_3</release>"),] #Contains other yeasts - maybe keep these? Probably not

for(i in 1:length(D2P2$Seqid)){
  D2P2$size[i]=D2P2$End[i]-D2P2$Start[i]    #Calculate the size of the disordered region
  print(i)
}
D2P2Processed=split(D2P2,D2P2$Seqid)    #Split the DF into a list where each entry is differentiated by the Seqid (the gene)
for(i in 1:length(D2P2Processed)){
  D2P2Processed[[i]]$Sequence[1] = as.character(SGD_Combined$Gene...Proteins...Sequence...Residues[which(SGD_Combined$Gene...Systematic.Name == as.character(D2P2Processed[[i]]$Seqid[[1]]))]) #Grab protein sequence for that gene from org.sc.sgd.db - this is useful for later on when we calculate prop disorder
  print(i)
}
dfD2P2Processed = (as.data.frame(do.call(rbind,D2P2Processed[])))
D2P2Proportions = lapply(D2P2Processed, function(x) aggregate(x$size,by=list(x$Predictor),function(i) c(totDisorder=sum(i), propDisorder=sum(i)/nchar(x$Sequence[1]), numDomains=sum(i>30),numDomainsStrict=sum(i>80), numDomainsSuper=sum(i>200))))#This takes each item in the list (in this case, each item is another list that contains a DF with info on that gene and the predicted regions) and sums up their disordered regions (totDisorder) and calculated the proportion of the protein that is predicted as disordered from that (propDisorder)
for(i in 1:length(D2P2Proportions)){
  D2P2Proportions[[i]]$Seqid=names(D2P2Proportions[i])  #Make a column with the seqID so we can group by that
}
dfD2P2Proportions=(as.data.frame(do.call(rbind,D2P2Proportions[]))) #Make it into a dataframe for easier plotting
dfD2P2Proportions=do.call(data.frame, dfD2P2Proportions)  #Remove matrix within dataframe



ggplot(dfD2P2Proportions, aes(x=Group.1, y = x.propDisorder))+geom_violin(alpha=0.2) +theme_Publication_GridlinesTesting()+scale_y_continuous(limits=c(0,1))
dfD2P2Proportions$ID="SGD"
bitr(dfD2P2Proportions$Seqid, fromType = "ORF", toType = c("GENENAME"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,3]
disorderNIP1 = dfD2P2Proportions[which(dfD2P2Proportions$Seqid %in% SGD2520$Gene...Systematic.Name),]
disorderNIP1$ID="NIP1"
disorderPDC1 = dfD2P2Proportions[which(dfD2P2Proportions$Seqid %in% SGD2521$Gene...Systematic.Name),]
disorderPDC1$ID="PDC1"
disorderCommon = dfD2P2Proportions[which(dfD2P2Proportions$Seqid %in% SGD_Common$Gene...Systematic.Name),]
disorderCommon$ID="Both"
disorderPropCombo=rbind(dfD2P2Proportions,disorderNIP1,disorderPDC1,disorderCommon)
disorderPropCombo$ID=factor(disorderPropCombo$ID,levels=c("SGD","NIP1","PDC1","Both"))

pDisorder = ggplot(disorderPropCombo,aes(x=ID,y=x.propDisorder))+
  geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_GridlinesTesting(base_family = "Avenir Black")+
  scale_y_continuous(limits=c(0,1))+ 
  facet_wrap(~disorderPropCombo$Group.1,ncol = 3)+ 
  scale_fill_manual(values=c("#939598","#FED208","#C697C4","#FF6F91")) +
  scale_colour_manual(values=c("#939598","#FEBB1A","#9C4993","#D63060")) +
  labs(y="Proportion of protein predicted to be disordered")

p1 <- ggplot_build(pDisorder)
#Now edit the plot data - the [[x]] accession is the order of plotting, so if you plot the beeswarm 1st or 3rd it would be [[1]] or [[3]], respectively
#We need to remove points that fall to the left of the centre line. Again, this depends on the number of x-axis factors you have.
#The below function takes points which fall to the left of the x-label (e.g. x-group is -ve) and makes these positive.
#It leaves values which have x-group as +ve because these don't need to be changed
p1$data[[1]]=p1$data[[1]] %>% mutate(x=case_when(
  x-round(x) < 0 ~ abs(x-round(x))+round(x),
  TRUE ~ x)
)

cairo_pdf("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder_Halfbox_Beeswarm.pdf",width=4.82,height=5.38)
plot_grid(ggplot_gtable(p1))
dev.off()

#Size IDR plot
ggplot(dfD2P2Processed, aes(x=log10(size)))+
  geom_histogram(aes(y=..density..), colour="black", fill="grey50",binwidth = 0.08)+
  geom_density(alpha=.2, fill=NA,col="red") +
  theme_Publication_Gridlines(base_family = "Avenir Black")+
  geom_vline(xintercept = log10(30), col="blue")+
  geom_vline(xintercept = log10(80), col = "purple")+
  geom_vline(xintercept = log10(200), col = "black")+
  facet_wrap(~Predictor)
  
#generate list of things with IDR >30aa and >80aa in more than 1 prediction method
#IDR >30bp and identified by more than one algorithm
IDRGenes = dfD2P2Proportions %>% 
  group_by(Seqid) %>% 
  summarise(IDRBinary=sum(x.numDomains>0),IDRBinaryStrict=sum(x.numDomainsStrict>0),IDRBinarySuper=sum(x.numDomainsSuper>0)) %>% 
  filter(IDRBinary>1) %>%
  data.frame()

#IDR >80bp and identified by more than one algorithm
IDRGenesStrict = dfD2P2Proportions %>% 
  group_by(Seqid) %>% 
  summarise(IDRBinary=sum(x.numDomains>0),IDRBinaryStrict=sum(x.numDomainsStrict>0),IDRBinarySuper=sum(x.numDomainsSuper>0)) %>% 
  filter(IDRBinaryStrict>1) %>%
  data.frame()

IDRGenesSuper = dfD2P2Proportions %>% 
  group_by(Seqid) %>% 
  summarise(IDRBinary=sum(x.numDomains>0),IDRBinaryStrict=sum(x.numDomainsStrict>0),IDRBinarySuper=sum(x.numDomainsSuper>0)) %>% 
  filter(IDRBinarySuper>1) %>%
  data.frame()

write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/IDRGenes.csv",IDRGenes)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/IDRGenesStrict.csv",IDRGenesStrict)
write.csv(file="~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/IDRGenesSuper.csv",IDRGenesSuper)

DisProt = read.table("~/Downloads/search_in_disprot.tsv",sep="\t",header = TRUE)
DisProtUnique = unique(subset(DisProt$acc, DisProt$confidence=="" & DisProt$end-DisProt$start > 80))
IDRUNIPROT = bitr(IDRGenesStrict$Seqid,fromType = "ENSEMBL", toType = "UNIPROT", OrgDb = org.Sc.sgd.db)[,2]

#New
IDRGenesStrict = dfD2P2Proportions %>% 
  group_by(Seqid, Start) %>% 
  summarise(IDRBinary=sum(x.numDomains>0),IDRBinaryStrict=sum(x.numDomainsStrict>0),IDRBinarySuper=sum(x.numDomainsSuper>0)) %>% 
  filter(IDRBinaryStrict>1) %>%
  data.frame()




#outliers = disorderPropCombo %>%
#  group_by(interaction(Group.1, ID)) %>%
#  mutate(outlier=ifelse(is_outlier(x.propDisorder) & ID != "SGD",CommonName, NA)) %>% filter(!is.na(outlier) & ID !="SGD") %>% summarise(outlier=count(outlier))

outliers = disorderPropCombo %>%
  group_by(Group.1, ID) %>%
  mutate(outlier=ifelse(x.propDisorder %in% unlist(list(boxplot.stats(x.propDisorder)$out)) & ID != "SGD",CommonName, NA)) %>% filter(!is.na(outlier)) %>% data.frame()

pDisorderLabelOutliers = disorderPropCombo %>%
  group_by(Group.1, ID) %>%
  mutate(outlier=ifelse(x.propDisorder %in% unlist(list(boxplot.stats(x.propDisorder)$out)) & ID != "SGD",CommonName, NA)) %>%
  ggplot(aes(x=ID,y=x.propDisorder))+
  geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  ggrepel::geom_text_repel(aes(label=outlier),family="Avenir Black", size=2.5)+
  geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_GridlinesTesting(base_family = "Avenir Black")+
  scale_y_continuous(limits=c(0,1))+ 
  facet_wrap(~disorderPropCombo$Group.1,ncol = 3)+ 
  scale_fill_manual(values=c("#939598","#FED208","#C697C4","#FF6F91")) +
  scale_colour_manual(values=c("#939598","#FEBB1A","#9C4993","#D63060")) +
  labs(y="Proportion of protein predicted to be disordered")


ggplot(disorderPropCombo,aes(x=ID,y=x.propDisorder))+geom_violin(alpha=0.2)+theme_Publication_GridlinesTesting()+scale_y_continuous(limits=c(0,1))+ facet_grid(disorderPropCombo$Group.1)
ggsave("~/SGA_BioInfo/Figs/Disorder.pdf",device = cairo_pdf,width=241, height=162, units="mm")
ggplot(disorderPropCombo,aes(x=ID,y=x.propDisorder))+geom_boxplot(aes(col=ID))+theme_Publication_GridlinesTesting()+scale_y_continuous(limits=c(0,1))+ facet_wrap(disorderPropCombo$Group.1,ncol = 3)+ scale_color_manual(values=c("#845EC2","#D65DB1","#FF6F91","#939598")) +labs(y="Proportion of protein predicted to be disordered")
ggsave("~/SGA_BioInfo/Figs/DisorderBoxplot.pdf",device = cairo_pdf,width = 160,height = 175,units = "mm")

#####PONDR#####
require(stringr)
require(gsubfn)
require(dplyr)
require(zoo)

readPondr <- function(pondr){
  pondr=readLines(pondr)
  delims = which(grepl(pattern = "==", pondr))
  delims = delims[c(1,2,3,5)]
  sampleName = gsub("=","",pondr[delims[1]])
  percDisordered = as.numeric(str_sub(pondr[(delims[2]+4)], -6,-1))
  disorderedRegions = data.frame(do.call(rbind, strapply(pondr[delims[2]+4:(delims[3]-3)], "\\[([^]]*)\\]", back = -1)))
  colnames(disorderedRegions) = c("Start","End")
  disorderedRegions$Perc = data.frame(Perc = do.call(rbind, lapply(pondr[(delims[2]+4):(delims[3]-1)], str_sub,start = -6,end=-1)))
  disorderBinary=paste(str_sub(pondr[delims[3]+seq(4,((delims[4]-delims[3])-6),4)], -55,-1),sep="",collapse="")
  disorderBinary=paste(disorderBinary,gsub(".*.\t "," ",str_sub(pondr[delims[4]-5],-54,-1)),sep="",collapse="")
  if(!is.integer(nchar(disorderBinary)/11)){
    short = 11-((nchar(disorderBinary)/11-as.integer(nchar(disorderBinary)/11))*11)
    fill=paste(rep(" ",short),collapse="")
    disorderBinary=paste(disorderBinary,fill,collapse="")
  }
  disorderBinary=paste(str_sub(disorderBinary, seq(from=2,to=(nchar(disorderBinary)-9),by=11), seq(from=11,to=nchar(disorderBinary),by=11)),collapse = "")
  detailedDisorder=data.frame(do.call(rbind,lapply(lapply(pondr[(delims[4]+3):length(pondr)], str_split, pattern="\t"), "[[",1)))
  colnames(detailedDisorder) = c("Residue","AminoAcid", "VSL2")
  for(i in 1:length(detailedDisorder[,1])){
    disorder = str_sub(disorderBinary,i,i)
    if(disorder == "D"){
      detailedDisorder$Binary[i] = 1
    }
    if(disorder==" "){
      detailedDisorder$Binary[i] = 0
    }
  }
  detailedDisorder$Residue = as.numeric(as.character(detailedDisorder$Residue))
  detailedDisorder$VSL2 = as.numeric(as.character(detailedDisorder$VSL2))
  return(list(sampleName,detailedDisorder,disorderedRegions))
}

plotPondr <-function(pondr, relativeX=FALSE, shareAxes=FALSE,k=75,ncol=1,strip.position="right"){
  if(!is.list(pondr) & relativeX==TRUE) {
    pondr[[2]]$Residue= (pondr[[2]]$Residue/max(pondr[[2]]$Residue))*100
  }
  if(is.list(pondr)&relativeX==TRUE){
    for(i in 1:length(pondr)){
      pondr[[i]][[2]]$Residue= (pondr[[i]][[2]]$Residue/max(pondr[[i]][[2]]$Residue))*100
    }
  }
  if(!is.list(pondr)&shareAxes==TRUE){
    pondrPlot = ggplot(pondr[[2]], aes(x=Residue, y=VSL2))+
      geom_line(col="grey50")+
      geom_line(data=pondr[[2]],aes(x=Residue, y=-Binary/20))+
      geom_area(data=pondr[[2]],aes(x=Residue, y=-Binary/20),fill="black")+
      theme_Publication_GridlinesTesting()+
      annotate("text", x=max(pondr[[2]]$Residue)-(max(pondr[[2]]$Residue)/30), y = 1, label=pondr[[1]], size = 3)+
      geom_line(data = data.frame(mean = rollmean(pondr[[2]]$VSL2,k = k)),aes(y=mean, x =pondr[[2]]$Residue[seq(1,length(pondr[[2]]$Residue),length(pondr[[2]]$Residue)/length(mean))]), col="red", size = 1.5)
    return(pondrPlot)
  }
  if(is.list(pondr)){
    pondr[[1]][[2]]$Name = pondr[[1]][[1]]
    bigDF = pondr[[1]][[2]]
    for(i in 2:length(pondr)){
      pondr[[i]][[2]]$Name = pondr[[i]][[1]]
      bigDF=rbind(bigDF, pondr[[i]][[2]])
    }
    if(shareAxes==TRUE){
      pondrPlot = ggplot(bigDF, aes(x=Residue, y=VSL2),group = Name)+
        geom_hline(yintercept = 0.5, col="cornflower blue", size = 0.75)+
        geom_line(col="grey50")+
        geom_line(data=bigDF,aes(x=Residue, y=-Binary/20))+
        geom_area(data=bigDF,aes(x=Residue, y=-Binary/20),fill="black")+
        theme_Publication_GridlinesTesting()+
        #annotate("text",data = aggregate(.~Name, bigDF, "mean"), x=max(bigDF$Residue)-(max(bigDF$Residue)/30), y = 1, label=bigDF$Name, group = Name,size = 3)+
        geom_line(data = cbind(do.call(rbind, lapply(split(bigDF[c(3,5)], bigDF$Name), function(x) data.frame(mean = rollmean(x[1][1], k=k)))),do.call(rbind,lapply(split(bigDF[c(1,2,4,5)], bigDF$Name), function(x) data.frame(x[seq(1,length(x$Name),(length(x$Name)/length(rollmean(x[1][1], k=k)))),])))),aes(y=VSL2, x =Residue, group=Name), col="red", size = 1.5)+
        facet_wrap(Name~.,ncol=ncol,strip.position = strip.position)
      return(pondrPlot)
    } 
    if(shareAxes==FALSE){
      pondrPlot = ggplot(bigDF, aes(x=Residue, y=VSL2),group = Name)+
        geom_hline(yintercept = 0.5, col="cornflower blue", size = 0.75)+
        geom_line(col="grey50")+
        geom_line(data=bigDF,aes(x=Residue, y=-Binary/20))+
        geom_area(data=bigDF,aes(x=Residue, y=-Binary/20),fill="black")+
        theme_Publication_GridlinesTesting()+
        #annotate("text",data = aggregate(.~Name, bigDF, "mean"), x=max(bigDF$Residue)-(max(bigDF$Residue)/30), y = 1, label=bigDF$Name, group = Name,size = 3)+
        geom_line(data = cbind(do.call(rbind, lapply(split(bigDF[c(3,5)], bigDF$Name), function(x) data.frame(mean = rollmean(x[1][1], k=k)))),do.call(rbind,lapply(split(bigDF[c(1,2,4,5)], bigDF$Name), function(x) data.frame(x[seq(1,length(x$Name),(length(x$Name)/length(rollmean(x[1][1], k=k)))),])))),aes(y=VSL2, x =Residue, group=Name), col="red", size = 1.5)+
        facet_wrap(Name~.,ncol=ncol,strip.position = strip.position, scales="free_x")
      return(pondrPlot)
    } 
  }
}

HTB1 = readPondr("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/HTB1_PondR.txt")
NCB2 = readPondr("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/NCB2_PondR.txt")
NHP2 = readPondr("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/NHP2_PondR.txt")
NPL3 = readPondr("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/NPL3_PondR.txt")
SLX8 = readPondr("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/SLX8_PondR.txt")
UME6 = readPondr("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/UME6_PondR.txt")
SPC19 = readPondr("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/SPC19_PondR.txt")
STB3 = readPondr("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/STB3_PondR.txt")
TIM11 = readPondr("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/TIM11_PondR.txt")

pondrList = list(HTB1, NCB2, NHP2, NPL3, SLX8, UME6, SPC19, STB3, TIM11)

cairo_pdf("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Disorder/Outliers_Disorder/plot.pdf", height = 5.83, width = 5.11)
plotPondr(pondrList,ncol = 3,strip.position = "top", k=40)
dev.off()

######Disorder stats######
#Using Man-whitney U test with bonferroni correction
#First need to split the dataframe into the diff algorithms
Espritz_N = subset(disorderPropCombo, disorderPropCombo$Group.1=="Espritz-N")
Espritz_D = subset(disorderPropCombo, disorderPropCombo$Group.1=="Espritz-D")
Espritz_X = subset(disorderPropCombo, disorderPropCombo$Group.1=="Espritz-X")
IUPred_L = subset(disorderPropCombo, disorderPropCombo$Group.1=="IUPred-L")
IUPred_S = subset(disorderPropCombo, disorderPropCombo$Group.1=="IUPred-S")
PrDOS = subset(disorderPropCombo, disorderPropCombo$Group.1=="PrDOS")
PV2 = subset(disorderPropCombo, disorderPropCombo$Group.1=="PV2")
VLXT = subset(disorderPropCombo, disorderPropCombo$Group.1=="VSL2b")

#Function to do test
disorderSigTest = function(df, mRNA){
  return(wilcox.test(x.propDisorder~ID, data=subset(df, df$ID=="SGD"|df$ID==mRNA)))
}

disorderSigTest(df=Espritz_N, mRNA = "NIP1")
disorderSigTest(df=Espritz_N, mRNA = "PDC1")
disorderSigTest(df=Espritz_N, mRNA = "Both")
disorderSigTest(df=Espritz_D, mRNA = "NIP1")
disorderSigTest(df=Espritz_D, mRNA = "PDC1")
disorderSigTest(df=Espritz_D, mRNA = "Both")
disorderSigTest(df=Espritz_X, mRNA = "NIP1")
disorderSigTest(df=Espritz_X, mRNA = "PDC1")
disorderSigTest(df=Espritz_X, mRNA = "Both")
disorderSigTest(df=IUPred_L, mRNA = "NIP1")
disorderSigTest(df=IUPred_L, mRNA = "PDC1")
disorderSigTest(df=IUPred_L, mRNA = "Both")
disorderSigTest(df=IUPred_S, mRNA = "NIP1")
disorderSigTest(df=IUPred_S, mRNA = "PDC1")
disorderSigTest(df=IUPred_S, mRNA = "Both")
disorderSigTest(df=PrDOS, mRNA = "NIP1")
disorderSigTest(df=PrDOS, mRNA = "PDC1")
disorderSigTest(df=PrDOS, mRNA = "Both")
disorderSigTest(df=PV2, mRNA = "NIP1")
disorderSigTest(df=PV2, mRNA = "PDC1")
disorderSigTest(df=PV2, mRNA = "Both")
disorderSigTest(df=VLXT, mRNA = "NIP1")
disorderSigTest(df=VLXT, mRNA = "PDC1")
disorderSigTest(df=VLXT, mRNA = "Both")


t#####PONDR#########
pondrList=list("~/SGA_BioInfo/Common_Disorder_Pondr/ARD1_PondR.txt", 
               "~/SGA_BioInfo/Common_Disorder_Pondr/BRP1_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/CHO2_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/GCD11_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/GPI17_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/HAL5_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/PHO2_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/PPT1_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/SAC1_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/SPC19_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/UME6_PondR.txt",
               "~/SGA_BioInfo/Common_Disorder_Pondr/YDJ1_PondR.txt")
pondrs=lapply(pondrList, "readPondr")
plotPondr(pondrs,k=10,ncol=3)
ggsave("~/SGA_BioInfo/Figs/Disorder_Common_Profiles_VSL2.pdf",device=cairo_pdf, width=241, height=162, units="mm")

########Plot the impacts of the deletions##########
catSGA=readxl::read_xlsx("~/SGA_BioInfo/Categorised_Hits_Pos_Neg_Misloc_CB.xlsx")
pPropCat = ggplot(catSGA) + 
  geom_bar(mapping=aes(x=Cat, y=..prop.., group = 1), stat="count")+
  facet_wrap(~mRNA)+
  theme_Publication_Gridlines()+
  scale_y_continuous(limits=c(0,1))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))

dev.off()
cairo_pdf("~/Documents/Thesis--sssss/Chapter2/4.1 A high-throughput SGA screen/ProportionCategoryHits.pdf", width = (108*0.0393701), height=(50*0.0393701))
pPropCat
dev.off()

######Comparisons with SG and PB data ######
PBSG = xlsx::read.xlsx("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/PBSG/Kershaw_PB_SG_Prots.xlsx",sheetIndex = 3)
PBSGReplete = xlsx::read.xlsx("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/PBSG/Kershaw_PB_SG_Prots.xlsx",sheetIndex = 4)
PB = PBSG[PBSG$GRANULE=="PB",]
SG = PBSG[PBSG$GRANULE=="SG",]
PBReplete = PBSGReplete[PBSGReplete$GRANULE=="PB",]
SGReplete = PBSGReplete[PBSGReplete$GRANULE=="SG",]

catPB = catSGA[catSGA$SeqID %in% PB$ORF,]
catSG = catSGA[catSGA$SeqID %in% SG$ORF,]
catPBReplete = catSGA[catSGA$SeqID %in% PBReplete$ORF,]
catSGReplete =catSGA[catSGA$SeqID %in% SGReplete$ORF,]

catSGA$PB = 0
catSGA[catSGA$SeqID %in% PB$ORF,]$PB=1
catSGA$SG = 0
catSGA[catSGA$SeqID %in% SG$ORF,]$SG=1
catSGA$PBReplete = 0
catSGA[catSGA$SeqID %in% PBReplete$ORF,]$PBReplete=1
catSGA$SGReplete = 0
catSGA[catSGA$SeqID %in% SGReplete$ORF,]$SGReplete=1



Reduce(intersect, list(catSGA[catSGA$mRNA == "NIP1",]$SeqID,PB$ORF,PBReplete$ORF))

draw.triple.venn(area1 = length(catSGA[catSGA$mRNA == "NIP1",]$SeqID), 
                 area2 = length(PB$ORF), 
                 area3 = length(PBReplete$ORF),
                 n12 = length(intersect(catSGA[catSGA$mRNA == "NIP1",]$SeqID,PB$ORF)),
                 n13 = length(intersect(catSGA[catSGA$mRNA == "NIP1",]$SeqID,PBReplete$ORF)),
                 n23 = length(intersect(PB$ORF,PBReplete$ORF)),
                 n123 = length(Reduce(intersect, list(catSGA[catSGA$mRNA == "NIP1",]$SeqID,PB$ORF,PBReplete$ORF))),
                 category=c("NIP1","PB","pre-PB"),
                 scaled=TRUE
)
draw.triple.venn(area1 = length(catSGA[catSGA$mRNA == "PDC1",]$SeqID), 
                 area2 = length(PB$ORF), 
                 area3 = length(PBReplete$ORF),
                 n12 = length(intersect(catSGA[catSGA$mRNA == "PDC1",]$SeqID,PB$ORF)),
                 n13 = length(intersect(catSGA[catSGA$mRNA == "PDC1",]$SeqID,PBReplete$ORF)),
                 n23 = length(intersect(PB$ORF,PBReplete$ORF)),
                 n123 = length(Reduce(intersect, list(catSGA[catSGA$mRNA == "PDC1",]$SeqID,PB$ORF,PBReplete$ORF))),
                 category=c("PDC1","PB","pre-PB"),
                 scaled=TRUE
)

draw.triple.venn(area1 = length(catSGA[catSGA$mRNA == "NIP1",]$SeqID), 
                 area2 = length(SG$ORF), 
                 area3 = length(SGReplete$ORF),
                 n12 = length(intersect(catSGA[catSGA$mRNA == "NIP1",]$SeqID,SG$ORF)),
                 n13 = length(intersect(catSGA[catSGA$mRNA == "NIP1",]$SeqID,SGReplete$ORF)),
                 n23 = length(intersect(SG$ORF,SGReplete$ORF)),
                 n123 = length(Reduce(intersect, list(catSGA[catSGA$mRNA == "NIP1",]$SeqID,SG$ORF,SGReplete$ORF))),
                 category=c("NIP1","SG","pre-SG")
)


draw.triple.venn(area1 = length(catSGA[catSGA$mRNA == "PDC1",]$SeqID), 
                 area2 = length(SG$ORF), 
                 area3 = length(SGReplete$ORF),
                 n12 = length(intersect(catSGA[catSGA$mRNA == "PDC1",]$SeqID,SG$ORF)),
                 n13 = length(intersect(catSGA[catSGA$mRNA == "PDC1",]$SeqID,SGReplete$ORF)),
                 n23 = length(intersect(SG$ORF,SGReplete$ORF)),
                 n123 = length(Reduce(intersect, list(catSGA[catSGA$mRNA == "PDC1",]$SeqID,SG$ORF,SGReplete$ORF))),
                 category=c("PDC1","SG","pre-SG"),
                 scaled=TRUE
)


####### Aggregated proteins in unstressed cells ######

agg = xlsx::read.xlsx("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Aggregates/Final analysis.xlsx", sheetIndex = 2)
agg$CommonName = bitr(agg$SGD, fromType = "ORF", toType = c("GENENAME"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,3]
agg2520 = subset(SGD2520, SGD2520$Gene...Standard.Name %in% agg$CommonName)
agg2521 = subset(SGD2521, SGD2521$Gene...Standard.Name %in% agg$CommonName)
aggInt = subset(SGD_Common, SGD2520$Gene...Standard.Name %in% agg$CommonName)
  
draw.pairwise.venn(area1 = length(SGD2520$Gene...Primary.DBID),
                   area2 = length(agg$SGD),
                   cross.area = length(intersect(SGD2520$Gene...Standard.Name, agg$CommonName)))
  
draw.pairwise.venn(area1 = length(SGD2521$Gene...Primary.DBID),
                   area2 = length(agg$SGD),
                   cross.area = length(intersect(SGD2521$Gene...Standard.Name, agg$CommonName)))

draw.triple.venn(area1 = length(SGD2520$Gene...Primary.DBID),
                 area2 = length(agg$SGD),
                 area3 = length(SGD2521$Gene...Primary.DBID),
                 n12 = length(intersect(SGD2520$Gene...Standard.Name, agg$CommonName)),
                 n23 = length(intersect(SGD2521$Gene...Standard.Name, agg$CommonName)),
                 n13 = length(intersect(SGD2520$Gene...Standard.Name, SGD2521$Gene...Standard.Name)),
                 n123 = length(Reduce(intersect, list(SGD2520$Gene...Standard.Name, SGD2521$Gene...Standard.Name,agg$CommonName))),
                 category = c("NIP1","Agg","PDC1"))
                   
                   
                   
