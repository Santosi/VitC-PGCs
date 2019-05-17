#####################################
#         PGC VitC RNAseq           #
#      Stephanie Parker 2017        #
#     sparkerditroia@gmial.com      #
#####################################

#bam files from Tophat alignment mm9 (-g 20)

library("Rsubread")
library("limma")
library("edgeR")
library("ggplot2")
library(VennDiagram)
library(BiocGenerics)
library(genefilter)
library(beadarray)
library(limma)
library(gdata)
library(ggplot2)
library(gplots)
library("Rsamtools")
library("DESeq2")
library(reshape2)
library(colorRamps)
library(rtracklayer)

setwd("~/Desktop/SPMP01 mm9 alignment/")

count_data <- featureCounts(files=c("1CF1_accepted_hits.bam","1CF2_accepted_hits.bam","2CF2_accepted_hits.bam","3CF2_accepted_hits.bam","4CF1_accepted_hits.bam","4CF2_accepted_hits.bam","1VF1_accepted_hits.bam","1VF2_accepted_hits.bam","2VF1_accepted_hits.bam","3VF1_accepted_hits.bam","4VF1_accepted_hits.bam","4VF2_accepted_hits.bam"),annot.ext="mp_genes.gtf",isGTFAnnotationFile=T,GTF.featureType="exon",GTF.attrType="gene_id",useMetaFeatures=T,allowMultiOverlap=F,strandSpecific=1,ignoreDup=T,reportReads=F,countMultiMappingReads=T,fraction=T)

count_data$counts[1:5,] #see first 1:5 rows

#-------------------------------------------------#
#      DESeq2 analysis
#-------------------------------------------------#
head(count_data)
counts <- count_data$counts
head(counts)
  #write.table(counts, "RawCounts_FeatureCounts_VitCRNAseq_Broad.txt",sep="\t", quote=F, row.names=T)
storage.mode(counts) = "integer"
round(colSums(counts)) #the total number of reads aligning to genes for each sample
head(cpm(counts["Stra8",])) #spot check known genes

#-----------------------------------------------------------------------------------------------------------------------------#
## Removal of genes with low expression
# only keep genes where >3 samples have a cpm >1. 
# Requires edgeR

is.expressed <- rowSums(cpm(counts)>1) >= 6
head(is.expressed)
  counts["Stra8",]
  counts[grep("Cyp26",rownames(counts)),]
counts_nocpm <- counts[-(is.expressed),]  #genes not expressed
head(counts_nocpm)
counts_cpm <- counts[is.expressed,]	#genes expressed in >3 samples
head(counts_cpm)

#-----------------------------------------------------------------------------------------------------------------------------#

coldata <- read.csv("ColData.csv",header = F, row.names=1)
condition <- factor(c("1CF1_accepted_hits.bam", "1CF2_accepted_hits.bam","2CF2_accepted_hits.bam",
                      "3CF2_accepted_hits.bam","4CF1_accepted_hits.bam","4CF2_accepted_hits.bam",
                      "1VF1_accepted_hits.bam","1VF2_accepted_hits.bam","2VF1_accepted_hits.bam",
                      "3VF1_accepted_hits.bam","4VF1_accepted_hits.bam","4VF2_accepted_hits.bam"))
matrix <- DESeqDataSetFromMatrix(countData=counts_cpm, colData=coldata, ~V2)
matrix
colnames(matrix)<-c("1CF1", "1CF2", "2CF2","3CF2","4CF1","4CF2","1VF1","1VF2","2VF1","3VF1","4VF1","4VF2")
exprs<-assay(matrix) #lets you look at the data (like exprs(X))
head(exprs)
summary(exprs)
round(colSums(exprs))

###########################  
### log transform and look
###########################

matrixR <- rlog(matrix)   #this is a variance stabilization method. it also normalises for sequencing depth, and log transforms.
exprsR<-assay(matrixR)   #ONLY USRE THIS FOR LOOKING AT, clustering, etc, USE RAW TABLE exprs FOR DE TESTING
head(exprsR)
summary(exprsR)
    write.table(exprsR, "VitaminC_RNAseq_norm-log2exprs.txt",sep="\t", quote=F, row.names=T)
exprsR<-as.data.frame(exprsR)


###########################
## check expression at select genes -- rlog variance stabilization of raw data (normalises for sequencing depth and log transforms)
###########################
exprsR["Dazl",] #check expression of something -- use to plot in prism / compare to qRT-PCR
exprsR["Sycp1",]
exprsR["Stra8",]
exprsR["Dnmt1",]
exprsR["Sars",] 

exprsR[grep("Cyp",rownames(exprsR)),]
exprsR[grep("Slc23",rownames(exprsR)),]
exprsR["Dppa3",]  # Stella / PGC7 - protect preimplantation imprints, expressed in oocytes, critical for post-fertilization 

dev.off() #reset graphics state
heatmap.2(as.matrix(exprsR[c("Smc1a","Msh2","Msh3","Msh4","Msh5","Msh6","Mlh3","Mlh1","Atm","Hormad1","Hormad2","Dmc1","Atr","Dazl","Sycp1","Sycp3","Sohlh2","Spo11","Taf7l","Stra8","Syce1","Mael","Gtsf1"),]), dendrogram="column",main="Expression of meiosis genes", col=colorpanel(40, "darkblue", "white","darkred"), key=TRUE, symkey=F, scale="row",density.info="none", trace="none",ColSideColors=colorBar)
# Rajkovic et al 2017, JCI, SOHLH1 and SOHLH2 coordinate oocyte differentiation


## Tet1(GT) and ESC+VitC genes
heatmap.2(as.matrix(exprsR[c("Dazl","Taf7l","Stra8","Hormad1","Mael","Sycp1","Sohlh2","Syce1"),]), dendrogram="column",main="VitaminC and Tet1 responsive meiosis genes", col=colorRampPalette(c("blue", "white","goldenrod1"))(n=299), key=TRUE, symkey=F, scale="row",density.info="none", trace="none",ColSideColors=colorBar,
          reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean))
#write.table(exprsR[c("Dazl","Taf7l","Stra8","Fkbp6","Hormad1","Mael","Sycp1","Sohlh2","Syce1"),],"Tet.ESC.genes.txt",sep="\t",quote=F, row.names=F)

## Tet expression and DMNT expression 
dev.off()
par(mar=c(7,4,4,2)+0.1) 
heatmap.2(as.matrix(exprsR[c("Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l","Tet1","Tet2","Tet3"),]),Colv=F,Rowv=F, dendrogram="none", col=colorRampPalette(c("blue", "white","goldenrod1"))(n=299), key=TRUE, symkey=F, scale="none",density.info="none", trace="none",ColSideColors=colorBar,
          reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),margins = c(10, 10),cexRow=2.5)
exprsR[c("Dnmt3a","Dnmt3b","Dnmt1","Tet1","Tet2","Tet3"),]
  write.table(exprsR[c("Dnmt3a","Dnmt3b","Dnmt1","Tet1","Tet2","Tet3"),],"Tet.DNMT.expresR.txt",sep="\t",quote=F, row.names=F)


###########################
## Correlation heatmaps 
###########################
pcorr<-function(x) as.dist(1-cor(t(x),method="pearson")) #make a separate variable for the distance matrix
data_corr<-cor(exprsR, method="pearson")  
data_corr
colorBar = c("azure3","azure3","azure3","azure3","azure3","azure3","brown3","brown3","brown3","brown3","brown3","brown3")
par(oma=c(3,4,3,4))
heatmap.2(data_corr, scale="none", col=colorpanel(20, "orange", "cornsilk"), key=TRUE, symkey=F, 
          density.info="none", trace="none",ColSideColors=colorBar, revC=T)

###########################
## Heatmap with 1000 most variable genes 
###########################
## uses package gene filter
##  rowVars = Row variance and standard deviation of a numeric array
head(exprsR)
topClusterGenes <- head( order( rowVars(as.matrix(exprsR )), decreasing=TRUE ), 1000) 
top1000GenesVar <- (exprsR[topClusterGenes,]) 
head(top1000GenesVar,10)
    # write.table(as.matrix(row.names(top500GenesVar)), "genes.txt",sep="\t", quote=F, row.names=F)
    # MOST VARIABLE GENE IS STRA8 AND 7TH IS SYCP1 
exprsRM<-as.matrix(exprsR)  
top_corr<-cor(exprsR[topClusterGenes,], method="pearson")  ##correlation with most variable genes
dev.off()
## correlation heatmap
heatmap.2(top_corr, scale="none", col=colorpanel(20, "orange", "cornsilk"), key=TRUE, symkey=F, 
          density.info="none", trace="none",ColSideColors=colorBar)
## gene expression heatmap
#my_palette2 <- colorRampPalette(c("goldenrod1", "cadetblue2", "cadetblue4"))(n=299)
my_palette3 <- colorRampPalette(c("cadetblue4", "white","goldenrod1"))(n=299)
dev.off()
heatmap.2(as.matrix(exprsR[topClusterGenes,]),scale="row", col=blue2red, cexRow=.1,
        trace="none",density.info="none",distfun=pcorr,ColSideColors=colorBar)
heatmap.2(as.matrix(exprsR[topClusterGenes,]),scale="row", col=colorRampPalette(c("blue", "white","goldenrod1"))(n=299), cexRow=.1,
          trace="none",density.info="none",distfun=pcorr,ColSideColors=colorBar)

###########################
## Only germ genes (identified by Michelle's PGC paper)
###########################
# correlations by germcell genes - according to MP data
germcellgenes_MP <- read.csv("GERMCELL_MP_TopTable-Female-padj05%2clogFC2.csv",header=T)
head(germcellgenes_MP)  # 8453 genes
germcelldata<-subset(exprsR,row.names(exprsR) %in% germcellgenes_MP[,1])  # 7836 in my data
topgerm <- head(order(rowVars(as.matrix(germcelldata)),decreasing = T),500)  # grab 500 most variable
top500germ <- (germcelldata[topgerm,])
germ_corr<-cor(germcelldata[topgerm,], method="pearson")
# heatmap on top 500 germ gene corr
par(oma=c(4,4,4,4))
heatmap.2(germ_corr, scale="none", col=colorpanel(20, "orange", "cornsilk"), key=TRUE, symkey=F, 
          density.info="none", trace="none")
# heatmap of gene expression
head(top500germ)
dev.off()
heatmap.2(as.matrix(top500germ),scale="row", col=my_palette3, cexRow=.1,
          main = "500 most variable germ genes", trace="none",density.info="none",distfun=pcorr)
#### TOP 50 germ ####
topgerm50 <- head(order(rowVars(as.matrix(germcelldata)),decreasing = T),50)  # grab 500 most variable
top50germ <- (germcelldata[topgerm50,])
germ_corr50<-cor(germcelldata[topgerm50,], method="pearson")
heatmap.2(germ_corr50, scale="none", col=colorpanel(20, "orange", "cornsilk"), key=TRUE, symkey=F, 
          density.info="none", trace="none",ColSideColors = colorBar)
head(top50germ)
dev.off()
heatmap.2(as.matrix(top50germ),scale="row", col=my_palette3, cexRow=.8,
          main = "50 most variable germ genes", trace="none",density.info="none",distfun=pcorr,ColSideColors = colorBar)

#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#


###########################
## DE expression changes -- DE analysis
###########################
head(matrix)
?DESeq
DEmatrix <- DESeq(matrix) #RAW data is processed to do Differential expression analysis NOT RLOG
head(DEmatrix)
exprsDE<-assay(DEmatrix)
head(exprsDE)

?results
toptable<-results(DEmatrix) #notes. P value is the Wald test p-value. 
mcols(toptable, use.names=T) #gives info on each column
plotMA(toptable, main="DESeq2", ylim=c(-3,3))
summary(toptable) 
      ##toptable_df <- as.data.frame(toptable)
      ##write.table(toptable_df, file="Pval_RNA.txt",sep="\t", quote=F, row.names=T)
toptableF<-subset(toptable, padj<0.1)
toptableF<-toptableF[order(toptableF$log2FoldChange),]
summary(toptableF)

## LOG FC > 1/-1 MEANS TWICE AS LARGE
toptableF_FC<-subset(toptableF, log2FoldChange < -1| log2FoldChange >1)    
toptableF_up<-subset(toptableF, log2FoldChange >.7 & padj < 0.05)    
toptableF_down<-subset(toptableF, log2FoldChange < -.7 & padj < 0.05)
upgenes<-as.data.frame(rownames(toptableF_up))      # 121 genes up -vitC p<.05 logFC<0.7 //  98 with 6 CMP>1
downgenes<-as.data.frame(rownames(toptableF_down))  # 390 genes down -vitC p<.05 logFC<-0.7 // 314 with 6 CMP>1
     write.table(upgenes, "up-98_6samplesCPM1.txt",sep="\t", quote=F, row.names=F)
     write.table(downgenes, "down-314_6samplesCPM1.txt",sep="\t", quote=F, row.names=F)

     
     ## up and down of MP germ genes
# germ_up<-subset(toptableF,row.names(toptableF) %in% germcellgenes_MP[,1]&log2FoldChange > .7)
#   upgenes.germ<-as.data.frame(rownames(germ_up))
# germ_down<-subset(toptableF,row.names(toptableF) %in% germcellgenes_MP[,1]&log2FoldChange< -0.7)
#   downgenes.germ<-as.data.frame(rownames(germ_down))

  
  ##############################
  ###### MA PLOT  #######
  ##############################
  
  toptable$Color = "grey14"
  toptable$Color[toptable$log2FoldChange > 0.7 & toptable$padj < 0.05] = "goldenrod"
  toptable$Color[toptable$log2FoldChange < -0.7 & toptable$padj< 0.05] = "blue"
  
  #toptable$Color[toptable$log2FoldChange > 0.7 & toptable$padj < 0.05] = "red"
  #toptable$Color[toptable$log2FoldChange < -0.7 & toptable$padj< 0.05] = "blue"
  #germ<-subset(toptable,row.names(toptable) %in% germcellgenes_MP[,1])
  #toptable$Color[row.names(toptable) %in% row.names(germ) & toptable$padj<0.05] = "grey"
  #imprint<-subset(toptable,row.names(toptable) %in% imprinted_genes[,1])
  #toptable$Color[row.names(toptable) %in% row.names(imprint)] = "purple"
  
  
  toptable2<-toptable[order(-(toptable$log2FoldChange)),]
  toptable2<-as.data.frame(toptable2)
  head(toptable2)
  
  sel<-toptable2[c("Dazl","Taf7l","Stra8","Fkbp6","Hormad1","Mael","Sycp1","Sohlh2","Syce1"),] #Tet1 germ genes
  sel
  sop<-toptable2[c("Chchd1","Sfrp1","Srp19","Rpusd1"),]  #most differential with targeted
  hk <- toptable2[c("Ubb","Rpl7","Fkbp6"),]  # qRT-PCR housekeeping genes
  ova <- toptable2[c("Xlr5b", "Xlr5a", "Tuba3b", "Tktl1", "Tcp11l2", "Taf7l", "Sycp1", "Stra8", "Slc25a31", "Nxf2", "Luzp4", "Gtsf1", "Ccdc73" ),]
            # genes up in VC ESC and down in Yamaguchi Tet1 KO
  imp <- toptable2[row.names(toptable) %in% row.names(imprint),]
    
    
  dev.off() #reset graphics state
  ggplot(data=toptable2, aes(log10(baseMean),log2FoldChange)) +       
    geom_point(alpha=.7, size=1, color=toptable2$Color) +
    geom_point(data=sel, aes(log10(baseMean),log2FoldChange), size=4, color="deeppink", stroke=1.2,shape=21) + 
    #geom_text(aes(data=sel,label=c("Dazl","Taf7l","Stra8","Fkbp6","Hormad1","Mael","Sycp1","Sohlh2","Syce1"),hjust=0, vjust=0)) +
    geom_hline(yintercept = 0, colour = "black", size=.2) +
    #geom_hline(yintercept = -1, colour = "black", size=1) +
    #geom_hline(yintercept = 1, colour = "red", size=0.8) +
    #geom_hline(yintercept = -1, colour = "blue", size=0.8) +
    geom_vline(xintercept = .8, colour = "black", size=.5) +
    ylim(c(-2, 2)) + xlim(c(0.8,6.5)) +
    xlab("log10 mean expression") + ylab("log2 fold change") +
    theme(panel.background=element_rect(fill="white"), panel.grid.major=element_line(color="white", size=1)) +
    theme(axis.text.x=element_text(size=20, color="black"), axis.text.y=element_text(size=20, color="black")) +
    ggtitle("Ctrl vs VitaminC removed E13.5 fPGCs - CPM > 1 6 samples") +
    theme(plot.title = element_text(size=30,lineheight=.8, vjust=2),axis.title.x=element_text(size=30),axis.title.y=element_text(size=30))
  

  
  
  ##############################
  ## Volcano plot
  ##############################
  plot(toptable$log2FoldChange, -log10(toptable$padj), col=toptable$Color, xlim = c(-2,2), 
       main="Volcano Plot", ylim = c(0,10), xlab="log2 fold change", ylab="-log10 pval", pch=20, cex=0.5)
  
  toptable$Color = "grey14"
  toptable$Color[toptable$log2FoldChange > .7 & toptable$padj < 0.05] = "red"
  toptable$Color[toptable$log2FoldChange < -.7 & toptable$padj< 0.05] = "blue"
  toptable3<-as.data.frame(toptable)
  head(toptable3)
  
  ggplot(data=toptable3, aes(log2FoldChange,-log10(padj)))+       
    geom_point(alpha=0.7, size=.8, color=toptable3$Color) +
    #geom_point(data=sel, aes(log2FoldChange,-log10(padj)), size=4, color="deeppink") +
    #geom_point(data=hk, aes(log2FoldChange,-log10(padj)), size=2, color="deeppink") +
    xlim(c(-2, 2))  +
    theme(panel.background=element_rect(fill="white"), panel.grid.major=element_line(color="white", size=1)) +
    theme(axis.text.x=element_text(size=20, color="black"), axis.text.y=element_text(size=20, color="black")) +
    #geom_hline(yintercept = 1.5, colour = "grey", size=0.3) +
    geom_vline(xintercept = .7, colour = "grey", size=0.3) +
    geom_vline(xintercept = -.7, colour = "grey", size=0.3) +
    ggtitle("Volcano") +
    theme(plot.title = element_text(size=30,lineheight=.8, vjust=2),axis.title.x=element_text(size=30),axis.title.y=element_text(size=30))
  
    
  
  
  
  #####################
  ###   PCA PLOT   ###
  #####################
  #USE PLOT PCA FROM THE DESEQ2 PACKAGE AND THEN USE GGPLOT2#
  PCA<-plotPCA(matrixR, intgroup = "V2", returnData=T, ntop=22000) #ntop=1000 top genes
  percentVar <- round(100 * attr(PCA, "percentVar")) 
  
  dev.off()
  qplot(PC1, PC2, color=group, data=PCA)+
    theme(text = element_text(size=36), axis.text.x = element_text(angle=90, vjust=1))+
    geom_point(size=4) +
    theme(legend.title=element_blank()) +
    theme(panel.grid.major=element_line(color="white", size=1)) +
    theme(axis.text.x=element_text(size=36, color="black"), axis.text.y=element_text(size=36, color="black")) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    scale_color_manual(values=c("azure4", "brown3"))+
    xlim(c(-25,25))+ylim(c(-25,25))
  
  
  
  



