############# MS List analysis #############
require(VennDiagram)
require(ggplot2)
require(dplyr)
require(clusterProfiler)
require(org.Sc.sgd.db)


theme_Publication_Gridlines <- function(axis_text_size=9, axis_title_size=12, base_family="Calibri",gridlines = FALSE,legendPos="bottom",legendDir = "horizontal") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background  = element_blank(),
            plot.background = element_rect(fill="transparent", colour=NA),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2,size = axis_title_size),
            axis.title.x = element_text(vjust = -0.2, size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
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




#Read in files#
ENO1MSList = read.csv("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO1.csv")
ENO2MSList = read.csv("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO2.csv")
NIP1MSList = read.csv("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/NIP1.csv")


ENO2BP = enrichGO(gene=ENO2MSList$GENENAME, OrgDb = org.Sc.sgd.db, ont = "BP", keyType = "GENENAME", pvalueCutoff = 0.05,universe = ENO2_Quantile_All[[2]]$GENENAME) #Enrichment analysis
ENO1BP = enrichGO(gene=ENO1MSList$GENENAME, OrgDb = org.Sc.sgd.db, ont = "BP", keyType = "GENENAME", pvalueCutoff = 0.05,universe = ENO1_Quantile_All[[2]]$GENENAME) #Enrichment analysis
NIP1BP = enrichGO(gene=NIP1MSList$GENENAME, OrgDb = org.Sc.sgd.db, ont = "BP", keyType = "GENENAME", pvalueCutoff = 0.05,universe = NIP1_Quantile_All[[2]]$GENENAME) #Enrichment analysis
ENO2CC = enrichGO(gene=ENO2MSList$GENENAME, OrgDb = org.Sc.sgd.db, ont = "CC", keyType = "GENENAME", pvalueCutoff = 0.05,universe = ENO2_Quantile_All[[2]]$GENENAME) #Enrichment analysis
ENO1CC = enrichGO(gene=ENO1MSList$GENENAME, OrgDb = org.Sc.sgd.db, ont = "CC", keyType = "GENENAME", pvalueCutoff = 0.05,universe = ENO1_Quantile_All[[2]]$GENENAME) #Enrichment analysis
NIP1CC = enrichGO(gene=NIP1MSList$GENENAME, OrgDb = org.Sc.sgd.db, ont = "CC", keyType = "GENENAME", pvalueCutoff = 0.05,universe = NIP1_Quantile_All[[2]]$GENENAME) #Enrichment analysis
ENO2MF = enrichGO(gene=ENO2MSList$GENENAME, OrgDb = org.Sc.sgd.db, ont = "MF", keyType = "GENENAME", pvalueCutoff = 0.05,universe = ENO2_Quantile_All[[2]]$GENENAME) #Enrichment analysis
ENO1MF = enrichGO(gene=ENO1MSList$GENENAME, OrgDb = org.Sc.sgd.db, ont = "MF", keyType = "GENENAME", pvalueCutoff = 0.05,universe = ENO1_Quantile_All[[2]]$GENENAME) #Enrichment analysis
NIP1MF = enrichGO(gene=NIP1MSList$GENENAME, OrgDb = org.Sc.sgd.db, ont = "MF", keyType = "GENENAME", pvalueCutoff = 0.05,universe = NIP1_Quantile_All[[2]]$GENENAME) #Enrichment analysis


#Venn Diagram#
draw.triple.venn(area1 = length(ENO1MSList$ORF),
                 area2 = length(ENO2MSList$ORF),
                 area3 = length(NIP1MSList$ORF),
                 n12 = length(intersect(ENO1MSList$ORF,ENO2MSList$ORF)),
                 n13 = length(intersect(ENO1MSList$ORF,NIP1MSList$ORF)),
                 n23 = length(intersect(ENO2MSList$ORF,NIP1MSList$ORF)),
                 n123 = length(Reduce(intersect, list(ENO1MSList$ORF,ENO2MSList$ORF,NIP1MSList$ORF))),
                 category=c("ENO1","ENO2","NIP1"),
                 euler.d=TRUE,
                 scaled = TRUE,
                 fontfamily = rep("Avenir Black",7),
                 cat.fontfamily=rep("Avenir Black", 3))

ENO1_int_ENO2 = intersect(ENO1MSList$ORF,ENO2MSList$ORF)
ENO1_int_NIP1 = intersect(ENO1MSList$ORF,NIP1MSList$ORF)
ENO2_int_NIP1 = intersect(ENO2MSList$ORF,NIP1MSList$ORF)
ENO1Only = Reduce(setdiff, list(ENO1MSList$ORF,ENO2MSList$ORF,NIP1MSList$ORF))
ENO2only = Reduce(setdiff, list(ENO2MSList$ORF,ENO1MSList$ORF,NIP1MSList$ORF))
NIP1only = Reduce(setdiff, list(NIP1MSList$ORF,ENO2MSList$ORF,ENO1MSList$ORF))

write.csv(file="~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO1_int_ENO2.csv",ENO1_int_ENO2)
write.csv(file="~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO1_int_NIP1.csv",ENO1_int_NIP1)
write.csv(file="~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO2_int_NIP1.csv",ENO2_int_NIP1)
write.csv(file="~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO1Only.csv",ENO1Only)
write.csv(file="~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO2only.csv",ENO2only)
write.csv(file="~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/NIP1Only.csv",NIP1only)

pantherENO1_int_ENO2 = read.table("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO1_int_ENO2PantherPie.txt",sep = "\t")
pantherENO1_int_NIP1 = read.table("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO1_int_NIP1PantherPie.txt",sep = "\t")
pantherNIP1_int_ENO2 = read.table("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO2_int_NIP1PantherPie.txt",sep = "\t")
pantherENO1 = read.table("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO1OnlyPantherPie.txt",sep = "\t")
pantherENO2 = read.table("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/ENO2OnlyPantherPie.txt",sep = "\t")
pantherNIP1 = read.table("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/QuantileNormalised/NIP1OnlyPantherPie.txt",sep = "\t")


pantherENO1_int_ENO2$ID = "ENO1IntersectENO2"
pantherENO1_int_NIP1$ID = "NIP1IntersectENO1"
pantherNIP1_int_ENO2$ID = "ENO2IntersectNIP1" 
pantherENO1$ID = "ENO1"
pantherENO2$ID = "ENO2"
pantherNIP1$ID = "NIP1"

pantherCombined = rbind(pantherENO1_int_ENO2,pantherENO1_int_NIP1,pantherNIP1_int_ENO2,pantherENO1,pantherENO2,pantherNIP1)

pantherCombined$ID = factor(as.factor(pantherCombined$ID),levels=c("ENO1","ENO2","NIP1","ENO1IntersectENO2","NIP1IntersectENO1","ENO2IntersectNIP1"))

write.csv(file="~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/PantherPieCollated.csv", pantherCombined)


pantherCombined$V5 = pantherCombined$V5 %>% substr(., start = 1,stop=nchar(as.character(.))-1) %>% as.numeric()
pantherCombined$V4 = pantherCombined$V4 %>% substr(., start = 1,stop=nchar(as.character(.))-1) %>% as.numeric()

PantherPlot = ggplot(pantherCombined, aes(x = ID, y = V5, fill = V2)) + 
  geom_bar(position = "fill",stat = "identity", col="black") +
  scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#FFC2E6", "#FFFFFF"))+
  theme_Publication_Gridlines(base_family = "Avenir Black",legendPos = "right",legendDir = "vertical")+
  labs(x="Dataset",y="Percentage of dataset",fill = "Panther class")

dev.off()
cairo_pdf("~/Documents/Thesis--sssss/Chapter3/blobs_clustering/October2020/PantherPlot.pdf", width = 8.06, height = 4.77)
PantherPlot
dev.off()




######Comparisons with SG and PB data ######
PBSG = xlsx::read.xlsx("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/PBSG/Kershaw_PB_SG_Prots.xlsx",sheetIndex = 3)
PBSGReplete = xlsx::read.xlsx("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/PBSG/Kershaw_PB_SG_Prots.xlsx",sheetIndex = 4)
PB = PBSG[PBSG$GRANULE=="PB",]
SG = PBSG[PBSG$GRANULE=="SG",]
PBReplete = PBSGReplete[PBSGReplete$GRANULE=="PB",]
SGReplete = PBSGReplete[PBSGReplete$GRANULE=="SG",]

ENO1MSList$mRNA = "ENO1"
ENO2MSList$mRNA = "ENO2"
NIP1MSList$mRNA = "NIP1"
FAPSList = rbind(ENO1MSList, ENO2MSList, NIP1MSList)

FAPSPB = FAPSList[FAPSList$ORF %in% PB$ORF,]
FAPSSG = FAPSList[FAPSList$ORF %in% SG$ORF,]
catPBReplete = FAPSList[FAPSList$ORF %in% PBReplete$ORF,]
catSGReplete =FAPSList[FAPSList$ORF %in% SGReplete$ORF,]

FAPSList$PB = 0
FAPSList[FAPSList$ORF %in% PB$ORF,]$PB=1
FAPSList$SG = 0
FAPSList[FAPSList$ORF %in% SG$ORF,]$SG=1
FAPSList$PBReplete = 0
FAPSList[FAPSList$ORF %in% PBReplete$ORF,]$PBReplete=1
FAPSList$SGReplete = 0
FAPSList[FAPSList$ORF %in% SGReplete$ORF,]$SGReplete=1



Reduce(intersect, list(FAPSList[FAPSList$mRNA == "NIP1",]$ORF,PB$ORF,PBReplete$ORF))

draw.triple.venn(area1 = length(FAPSList[FAPSList$mRNA == "NIP1",]$ORF), 
                 area2 = length(PB$ORF), 
                 area3 = length(PBReplete$ORF),
                 n12 = length(intersect(FAPSList[FAPSList$mRNA == "NIP1",]$ORF,PB$ORF)),
                 n13 = length(intersect(FAPSList[FAPSList$mRNA == "NIP1",]$ORF,PBReplete$ORF)),
                 n23 = length(intersect(PB$ORF,PBReplete$ORF)),
                 n123 = length(Reduce(intersect, list(FAPSList[FAPSList$mRNA == "NIP1",]$ORF,PB$ORF,PBReplete$ORF))),
                 category=c("NIP1","PB","pre-PB"),
                 fontfamily = rep("Avenir Black",7),
                 cat.fontfamily=rep("Avenir Black", 3),
                 scaled=TRUE
)
draw.triple.venn(area1 = length(FAPSList[FAPSList$mRNA == "ENO1",]$ORF), 
                 area2 = length(PB$ORF), 
                 area3 = length(PBReplete$ORF),
                 n12 = length(intersect(FAPSList[FAPSList$mRNA == "ENO1",]$ORF,PB$ORF)),
                 n13 = length(intersect(FAPSList[FAPSList$mRNA == "ENO1",]$ORF,PBReplete$ORF)),
                 n23 = length(intersect(PB$ORF,PBReplete$ORF)),
                 n123 = length(Reduce(intersect, list(FAPSList[FAPSList$mRNA == "ENO1",]$ORF,PB$ORF,PBReplete$ORF))),
                 category=c("ENO1","PB","pre-PB"),
                 fontfamily = rep("Avenir Black",7),
                 cat.fontfamily=rep("Avenir Black", 3),
                 scaled=TRUE
)

draw.triple.venn(area1 = length(FAPSList[FAPSList$mRNA == "ENO2",]$ORF), 
                 area2 = length(PB$ORF), 
                 area3 = length(PBReplete$ORF),
                 n12 = length(intersect(FAPSList[FAPSList$mRNA == "ENO2",]$ORF,PB$ORF)),
                 n13 = length(intersect(FAPSList[FAPSList$mRNA == "ENO2",]$ORF,PBReplete$ORF)),
                 n23 = length(intersect(PB$ORF,PBReplete$ORF)),
                 n123 = length(Reduce(intersect, list(FAPSList[FAPSList$mRNA == "ENO2",]$ORF,PB$ORF,PBReplete$ORF))),
                 category=c("ENO2","PB","pre-PB"),
                 fontfamily = rep("Avenir Black",7),
                 cat.fontfamily=rep("Avenir Black", 3),
                 scaled=TRUE
)

draw.triple.venn(area1 = length(FAPSList[FAPSList$mRNA == "NIP1",]$ORF), 
                 area2 = length(SG$ORF), 
                 area3 = length(SGReplete$ORF),
                 n12 = length(intersect(FAPSList[FAPSList$mRNA == "NIP1",]$ORF,SG$ORF)),
                 n13 = length(intersect(FAPSList[FAPSList$mRNA == "NIP1",]$ORF,SGReplete$ORF)),
                 n23 = length(intersect(SG$ORF,SGReplete$ORF)),
                 n123 = length(Reduce(intersect, list(FAPSList[FAPSList$mRNA == "NIP1",]$ORF,SG$ORF,SGReplete$ORF))),
                 fontfamily = rep("Avenir Black",7),
                 cat.fontfamily=rep("Avenir Black", 3),
                 category=c("NIP1","SG","pre-SG")
)


draw.triple.venn(area1 = length(FAPSList[FAPSList$mRNA == "ENO1",]$ORF), 
                 area2 = length(SG$ORF), 
                 area3 = length(SGReplete$ORF),
                 n12 = length(intersect(FAPSList[FAPSList$mRNA == "ENO1",]$ORF,SG$ORF)),
                 n13 = length(intersect(FAPSList[FAPSList$mRNA == "ENO1",]$ORF,SGReplete$ORF)),
                 n23 = length(intersect(SG$ORF,SGReplete$ORF)),
                 n123 = length(Reduce(intersect, list(FAPSList[FAPSList$mRNA == "ENO1",]$ORF,SG$ORF,SGReplete$ORF))),
                 fontfamily = rep("Avenir Black",7),
                 cat.fontfamily=rep("Avenir Black", 3),
                 category=c("ENO1","SG","pre-SG"),
                 scaled=TRUE
)

draw.triple.venn(area1 = length(FAPSList[FAPSList$mRNA == "ENO2",]$ORF), 
                 area2 = length(SG$ORF), 
                 area3 = length(SGReplete$ORF),
                 n12 = length(intersect(FAPSList[FAPSList$mRNA == "ENO2",]$ORF,SG$ORF)),
                 n13 = length(intersect(FAPSList[FAPSList$mRNA == "ENO2",]$ORF,SGReplete$ORF)),
                 n23 = length(intersect(SG$ORF,SGReplete$ORF)),
                 n123 = length(Reduce(intersect, list(FAPSList[FAPSList$mRNA == "ENO2",]$ORF,SG$ORF,SGReplete$ORF))),
                 category=c("ENO2","SG","pre-SG"),
                 fontfamily = rep("Avenir Black",7),
                 cat.fontfamily=rep("Avenir Black", 3)
)


####### Aggregated proteins in unstressed cells ######

agg = xlsx::read.xlsx("~/Documents/Thesis--sssss/Chapter2/4.2 Analysis of protein hits/Aggregates/Final analysis.xlsx", sheetIndex = 2)
agg$CommonName = bitr(agg$SGD, fromType = "ORF", toType = c("GENENAME"), OrgDb = org.Sc.sgd.db,drop = FALSE)[,3]
aggFAPS = subset(FAPSList, FAPSList$ORF %in% agg$SGD)

draw.quad.venn(area1 = length(ENO1MSList$ORF),
               area2 = length(ENO2MSList$ORF),
               area3 = length(NIP1MSList$ORF),
               area4 = length(agg$SGD),
               n12 = length(intersect(ENO1MSList$ORF,ENO2MSList$ORF)),
               n13 = length(intersect(ENO1MSList$ORF,NIP1MSList$ORF)),
               n14 = length(intersect(ENO1MSList$ORF,agg$SGD)),
               n23 = length(intersect(ENO2MSList$ORF,NIP1MSList$ORF)),
               n24 = length(intersect(ENO2MSList$ORF,agg$SGD)),
               n34 = length(intersect(NIP1MSList$ORF,agg$SGD)),
               n123 = length(Reduce(intersect, list(ENO1MSList$ORF,ENO2MSList$ORF,NIP1MSList$ORF))),
               n124 = length(Reduce(intersect, list(ENO1MSList$ORF,ENO2MSList$ORF,agg$SGD))),
               n134 = length(Reduce(intersect, list(ENO1MSList$ORF,NIP1MSList$ORF,agg$SGD))),
               n234 = length(Reduce(intersect, list(ENO2MSList$ORF,NIP1MSList$ORF,agg$SGD))),
               n1234 = length(Reduce(intersect, list(ENO1MSList$ORF,ENO2MSList$ORF,NIP1MSList$ORF,agg$SGD))),
               category=c("ENO1","ENO2","NIP1","Aggregates"),
               fontfamily = rep("Avenir Black",15),
               cat.fontfamily=rep("Avenir Black", 4))
  
  
#Fischer's exact
fisher.test(matrix(c(length(intersect(NIP1MSList$ORF,agg$SGD)),length(setdiff(agg$SGD,NIP1MSList$ORF)),length(NIP1MSList$ORF),length(agg$SGD)),nrow=2,ncol=2),alternative="greater")

m = length(agg$SGD)
N=length(NIP1_Quantile_All[[2]]$ORF) #Total number of genes in yeast (in an ideal world, this number would be less, e.g. the number of proteins potentially discoverable by the methods)
n=N-m
k = length(NIP1MSList$ORF)
x = length(intersect(NIP1MSList$ORF,agg$SGD))
p.valueNIP1=phyper(q=x-1,m=m,n=n,k=k,lower.tail = FALSE)

m = length(agg$SGD)
N=length(ENO1_Quantile_All[[2]]$ORF) #Total number of genes in yeast (in an ideal world, this number would be less, e.g. the number of proteins potentially discoverable by the methods)
n=N-m
k = length(ENO1MSList$ORF)
x = length(intersect(ENO1MSList$ORF,agg$SGD))
p.valueENO1=phyper(q=x-1,m=m,n=n,k=k,lower.tail = FALSE)

m = length(agg$SGD)
N=length(ENO2_Quantile_All[[2]]$ORF) #Total number of genes in yeast (in an ideal world, this number would be less, e.g. the number of proteins potentially discoverable by the methods)
n=N-m
k = length(ENO2MSList$ORF)
x = length(intersect(ENO2MSList$ORF,agg$SGD))
p.valueENO2=phyper(q=x-1,m=m,n=n,k=k,lower.tail = FALSE)

write.csv(file="~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/Aggregates/StatTestEnrichment.csv", data.frame(Condition=c("NIP1","ENO1","ENO2"), pVal = c(p.valueNIP1, p.valueENO1,p.valueENO2)))


######Making FASTA for slider
ENO1Fasta = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/ENO1_FASTA.csv", header = FALSE)
ENO1ALLFasta = read.table("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/ENO1All_FASTA.tsv", sep="\t", header=FALSE)
ENO2Fasta = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/ENO2_FASTA.csv",header=FALSE)
ENO2ALLFasta = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/ENO2All_FASTA.csv",header=FALSE)
NIP1Fasta = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/NIP1_FASTA.csv",header=FALSE)
NIP1ALLFasta = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/NIP1All_FASTA.csv",header=FALSE)

write.table(file="~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/ENO1_FASTA.txt", data.frame(paste(">",ENO1Fasta$V4,sep=""),stringr::str_replace(ENO1Fasta$V6, pattern = "\\*","")),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\n")
write.table(file="~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/ENO1ALL_FASTA.txt", data.frame(paste(">",ENO1ALLFasta$V4,sep=""),stringr::str_replace(ENO1ALLFasta$V6, pattern = "\\*","")),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\n")
write.table(file="~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/ENO2_FASTA.txt", data.frame(paste(">",ENO2Fasta$V4,sep=""),stringr::str_replace(ENO2Fasta$V6, pattern = "\\*","")),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\n")
write.table(file="~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/ENO2ALL_FASTA.txt", data.frame(paste(">",ENO2ALLFasta$V4,sep=""),stringr::str_replace(ENO2ALLFasta$V6, pattern = "\\*","")),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\n")
write.table(file="~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/NIP1_FASTA.txt", data.frame(paste(">",NIP1Fasta$V4,sep=""),stringr::str_replace(NIP1Fasta$V6, pattern = "\\*","")),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\n")
write.table(file="~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/NIP1ALL_FASTA.txt", data.frame(paste(">",NIP1ALLFasta$V4,sep=""),stringr::str_replace(NIP1ALLFasta$V6, pattern = "\\*","")),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\n")

####SLIDER#####
ENO1Slider = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/ENO1SLIDER.csv")
ENO1ALLSlider = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/ENO1ALLSLIDER.csv")
ENO2Slider = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/ENO2SLIDER.csv")
ENO2ALLSlider = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/ENO2ALLSLIDER.csv")
NIP1Slider = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/NIP1SLIDER.csv")
NIP1ALLSlider = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/NIP1ALLSLIDER.csv")
ProteomeSlider = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/CerevisiaeProteome_SLIDER.csv")
mRNPSlider = read.csv("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/SliderResults_Common.csv")

ENO1Slider$ID = "ENO1"
ENO1ALLSlider$ID = "ENO1ALL"
ENO2Slider$ID = "ENO2"
ENO2ALLSlider$ID = "ENO2ALL"
NIP1Slider$ID = "NIP1"
NIP1ALLSlider$ID ="NIP1ALL"
ProteomeSlider$ID ="BKGD"
mRNPSlider$ID = "IntmRNP"

SliderCombo = rbind(ENO1Slider, ENO1ALLSlider, ENO2Slider, ENO2ALLSlider, NIP1Slider, NIP1ALLSlider, ProteomeSlider, mRNPSlider)

SliderCombo$ID =   factor(as.factor(SliderCombo$ID),levels=c("BKGD","ENO1ALL","ENO1","ENO2ALL","ENO2","NIP1ALL","NIP1","IntmRNP"))

pDisorder = ggplot(SliderCombo, aes(x=ID, y=SLIDER.score..the.higher.the.more.likely.a.protein.has.long.disorder.segment.))+
  ggbeeswarm::geom_quasirandom(aes(col=ID,fill=ID),shape=21)+
  gghalves::geom_half_boxplot(aes(col=ID), fill="white")+
  theme_Publication_Gridlines(base_family = "Avenir Black")+
  scale_fill_manual(values=c("#C2C3C6","#F4B368","#D88D42","#F4B368","#D88D42","#F16058","#B84845","#3F59A8")) +
  scale_colour_manual(values=c("#C2C3C6","#F4B368","#D88D42","#F4B368","#D88D42","#F16058","#B84845","#3F59A8")) +
  labs(y="Slider Disorder Score", x="Dataset")+
  theme(legend.position="none")+
  theme(plot.margin = unit(c(1, 3, 1, 1), "mm"))

p1 <- ggplot_build(pDisorder)
#Now edit the plot data - the [[x]] accession is the order of plotting, so if you plot the beeswarm 1st or 3rd it would be [[1]] or [[3]], respectively
#We need to remove points that fall to the left of the centre line. Again, this depends on the number of x-axis factors you have.
#The below function takes points which fall to the left of the x-label (e.g. x-group is -ve) and makes these positive.
#It leaves values which have x-group as +ve because these don't need to be changed
p1$data[[1]]=p1$data[[1]] %>% mutate(x=case_when(
  x-round(x) < 0 ~ abs(x-round(x))+round(x),
  TRUE ~ x)
)

cairo_pdf("~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/SliderPlot.pdf", width = (120*0.039),height = (80*0.039))
cowplot::plot_grid(ggplot_gtable(p1))
dev.off()

test1=wilcox.test(SLIDER.score..the.higher.the.more.likely.a.protein.has.long.disorder.segment.~ID, data=subset(SliderCombo, ID=="BKGD"|ID=="ENO1ALL"))
test2=wilcox.test(SLIDER.score..the.higher.the.more.likely.a.protein.has.long.disorder.segment.~ID, data=subset(SliderCombo, ID=="ENO1ALL"|ID=="ENO1"))
test3=wilcox.test(SLIDER.score..the.higher.the.more.likely.a.protein.has.long.disorder.segment.~ID, data=subset(SliderCombo, ID=="BKGD"|ID=="ENO2ALL"))
test4=wilcox.test(SLIDER.score..the.higher.the.more.likely.a.protein.has.long.disorder.segment.~ID, data=subset(SliderCombo, ID=="ENO2ALL"|ID=="ENO2"))
test5=wilcox.test(SLIDER.score..the.higher.the.more.likely.a.protein.has.long.disorder.segment.~ID, data=subset(SliderCombo, ID=="BKGD"|ID=="NIP1ALL"))
test6=wilcox.test(SLIDER.score..the.higher.the.more.likely.a.protein.has.long.disorder.segment.~ID, data=subset(SliderCombo, ID=="NIP1ALL"|ID=="NIP1"))
test7=wilcox.test(SLIDER.score..the.higher.the.more.likely.a.protein.has.long.disorder.segment.~ID, data=subset(SliderCombo, ID=="BKGD"|ID=="IntmRNP"))

write.csv(file="~/Documents/Thesis--sssss/Chapter3/Figure5.5_FAPS-MS-ProteinAnalysis/PB_SG/SLIDER/Stats.csv", data.frame(comp=c("ENO1ALL","ENO1","ENO2ALL","ENO2","NIP1ALL","NIP1","IntmRNP"), pVal = c(test1$p.value, test2$p.value, test3$p.value, test4$p.value, test5$p.value, test6$p.value, test7$p.value)))


########Co-trans chaperones########
SSB = xlsx::read.xlsx("~/Documents/Thesis--sssss/Chapter3/Figure5.6_CotransFolding/Ssb1_2_Targets.xlsx", sheetIndex = 1)
NIP1SSB = read.csv
localizedRNAs = localizedRNAs %>% 
  mutate(SSBTarget = case_when(ORF %in% SSB$systematic.name ~ 1,
  TRUE ~ 0)
  )
SSBlocalised = localizedRNAs %>% group_by(Granule) %>% summarise(prop=sum(SSBTarget)/n())
SSBlocalised=rbind(data.frame(Granule="BKGD", prop=length(unique(SSB$systematic.name))/6275),SSBlocalised)

SSBlocalised$Granule =  factor(as.factor(SSBlocalised$Granule ),levels=c("BKGD","CoFe","Afe","Translation"))

cairo_pdf("~/Documents/Thesis--sssss/Chapter3/Figure5.6_CotransFolding/SSB1_2_EnrichmentHyperGeom.pdf", width = (90*0.039), height = (90*0.039))
ggplot(SSBlocalised, aes(x=Granule, y=prop))+
  geom_bar(stat = "identity")+
  theme_Publication_Gridlines(base_family = "Avenir Black")
dev.off()

#Enrichment
m=length(unique(SSB$systematic.name)) #Number of marks, e.g. the total number of genes with RBP function annotated
N=6275 #Total number of genes in yeast (in an ideal world, this number would be less, e.g. the number of proteins potentially discoverable by the methods)
n=N-m #non-marked elements
k=length(subset(localizedRNAs$mRNA, localizedRNAs$Granule=="CoFe")) #Background number of genes we were analysing - e.g. the number of genes identified in 2520
x=length(subset(localizedRNAs$mRNA, localizedRNAs$Granule=="CoFe" & localizedRNAs$SSBTarget==1)) #Number of marked genes in background
p.valueCoFe=phyper(q=x-1,m=m,n=n,k=k,lower.tail = FALSE)

m=length(unique(SSB$systematic.name)) #Number of marks, e.g. the total number of genes with RBP function annotated
N=6275 #Total number of genes in yeast (in an ideal world, this number would be less, e.g. the number of proteins potentially discoverable by the methods)
n=N-m #non-marked elements
k=length(subset(localizedRNAs$mRNA, localizedRNAs$Granule=="Afe")) #Background number of genes we were analysing - e.g. the number of genes identified in 2520
x=length(subset(localizedRNAs$mRNA, localizedRNAs$Granule=="Afe" & localizedRNAs$SSBTarget==1)) #Number of marked genes in background
p.valueAFe=phyper(q=x-1,m=m,n=n,k=k,lower.tail = FALSE)

m=length(unique(SSB$systematic.name)) #Number of marks, e.g. the total number of genes with RBP function annotated
N=6275 #Total number of genes in yeast (in an ideal world, this number would be less, e.g. the number of proteins potentially discoverable by the methods)
n=N-m #non-marked elements
k=length(subset(localizedRNAs$mRNA, localizedRNAs$Granule=="Translation")) #Background number of genes we were analysing - e.g. the number of genes identified in 2520
x=length(subset(localizedRNAs$mRNA, localizedRNAs$Granule=="Translation" & localizedRNAs$SSBTarget==1)) #Number of marked genes in background
p.valueTranslation=phyper(q=x-1,m=m,n=n,k=k,lower.tail = FALSE)
  
write.csv(file="~/Documents/Thesis--sssss/Chapter3/Figure5.6_CotransFolding/SSB1_2_EnrichmentHyperGeom.csv", data.frame(Condition=c("CoFe","AFe","Translation"), pVal = c(p.valueCoFe, p.valueAFe,p.valueTranslation)))


