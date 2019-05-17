############################################
######        E13.5 PGC RRBS          ######
######    DMR detection with BiSeq    ######
######     CpGs with >5x coverage     ######
######         in 9/12 samples        ######
######      Stephanie Parker 2016     ######

## BiSeq reference: 
## Hebestreit, K., Dugas, M. & Klein, H.-U. 
## Detection of significantly differentially 647 methylated regions in 
## targeted bisulfite sequencing data. Bioinformatics 29, 1647–1653 648 (2013). 649

library(BiSeq)  
library(gplots)
library(ggplot2)
library(colorRamps)
library(dichromat)
library(genefilter)
library(rtracklayer)
library(reshape2)
library(BiSeq)  
library(colorRamps)

setwd("~/Desktop/RRBSAnalysis_SLP01_mm9Broad/9:12 samples/")

#-----------------------------------------------------------------------------------------------------------------------------#

## List files and create BSraw object
files <- list.files(pattern = "_mm9_rrbs.bedgraph.bismark.cov")
files
bisseq <- readBismark(files, 
                      colData = DataFrame(group=c("ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","-vitc","-vitc","-vitc","-vitc","-vitc","-vitc"),
                                      row.names=c("C1F3","C2F2","C2F3","C3F1","C3F2","C3F3","V1F2","V2F1","V2F2","V3F1","V3F2","V4F1")))
covStatistics(bisseq)              
covBoxplots(bisseq,ylim=c(0,100))


## Filter by coverage -- 5x across all samples
bisseq.cov <- filterByCov(bisseq,minCov=5,global=F)  
covStatistics(bisseq.cov)
head(totalReads(bisseq.cov))


## 90th quantile of coverage to remove high coverage bias
bs.cov <- totalReads(bisseq.cov) > 0
quantile(totalReads(bisseq.cov)[bs.cov], 0.9)
## 27x coverage is the cutoff for 90th percentile in this cohort ##
covBoxplots(bisseq.cov,ylim=c(0,100),main="90th percentile coverage cutoff - 27x")
abline(h=23,col="red")
cov.lim <- limitCov(bisseq.cov, maxCov = 27) 
covBoxplots(cov.lim)
covStatistics(cov.lim)    
	# 590146 CpGs covered 5x in all samples
raw.methylation <- rawToRel(cov.lim)

#-----------------------------------------------------------------------------------------------------------------------------#

### Viewing data ###
colData(cov.lim)
head(totalReads(cov.lim))
head(methReads(cov.lim))
rrbs.rel <- rawToRel(cov.lim)
head(methLevel(rrbs.rel))
barplot(colMeans(methLevel(rrbs.rel),na.rm=T),col=colorBar, ylab = "methylation")

#-----------------------------------------------------------------------------------------------------------------------------#

### Correlation plots ### 
corr<-cor(methLevel(rrbs.rel),use = "na.or.complete", method="pearson")
heatmap.2(corr, scale="none", col=colorpanel(20, "orange", "cornsilk"), key=TRUE, symkey=F, na.rm=T,
          density.info="none", trace="none",ColSideColors=colorBar)

#-----------------------------------------------------------------------------------------------------------------------------#


# Subset DAZL or SYCP1 #
## Define genome region and subset methylation data
colorBar = c("azure3","azure3","azure3","azure3","azure3","azure3","brown3","brown3","brown3","brown3","brown3","brown3")
# dazl.region.interest <- GRanges(seqnames = "chr17",
#                                 ranges=IRanges(start = 50432598, 
#                                                end = 50432846))
# Below is Dazl promoter only
dazl.region.interest <- GRanges(seqnames = "chr17",
                                ranges=IRanges(start = 50431059, 
                                               end = 50434856))
sub.dazl.interest <- subsetByOverlaps(cov.lim, dazl.region.interest)
totalReads(sub.dazl.interest)
methReads(sub.dazl.interest)
covStatistics(sub.dazl.interest)
sub.dazl.interest.rel <- rawToRel(sub.dazl.interest)
dazl<-methLevel(sub.dazl.interest.rel)
dazl
boxplot(methLevel(sub.dazl.interest.rel),col=colorBar,main="Dazl")
# no dendrogram
?dichromat
?colorschemes
heatmap.2(dazl,na.color="black",scale="none", col=colorschemes$DarkRedtoBlue.18, cexRow=.1,
          trace="none",dendrogram='none', Rowv=FALSE, Colv=FALSE,density.info="none",ColSideColors=colorBar)
# cluster
heatmap.2(dazl,na.color="black",scale="none", col=colorschemes$DarkRedtoBlue.18, cexRow=.1,
          trace="none",dendrogram="column",Rowv=FALSE, Colv=T,density.info="none",ColSideColors=colorBar)


sycp1.region <- GRanges(seqnames = "chr3",
                        ranges=IRanges(start = 102729640, 
                                       end = 102749950))
sub.sycp1 <- subsetByOverlaps(cov.lim, sycp1.region)
totalReads(sub.sycp1)
methReads(sub.sycp1)
covStatistics(sub.sycp1)
sub.sycp1.rel <- rawToRel(sub.sycp1)
sycp1<-methLevel(sub.sycp1.rel)
boxplot(methLevel(sub.sycp1.rel),col=colorBar, main="Sycp1")
sycp1 <- sycp1[rowVars(sycp1,na.rm=T)!=0,] # remove rows with zero variability 
dev.off()
heatmap.2(sycp1,scale="none", col=green2red, cexRow=.1,na.rm=T,dendrogram="column",Rowv=FALSE,
          trace="none",density.info="none",ColSideColors=colorBar)
heatmap.2(sycp1,scale="row", col=colorschemes$DarkRedtoBlue.18, cexRow=.1,na.rm=T,dendrogram="column",Rowv=FALSE,
          trace="none",density.info="none",ColSideColors=colorBar)

#-----------------------------------------------------------------------------------------------------------------------------#

## GENERATING CLUSTERS & SMOOTHING ##

## Cluster CpGs that are covered >5x (could use 9/12 if bisseq used)
## Minimum CpGs per cluster = 5
## Minimum distance between CpGs = 20bp
clusters <- clusterSites(object = cov.lim,
                               perc.samples = 9/12,
                               min.sites = 5,
                               max.dist = 20)
covStatistics(clusters)
granges(clusters)
# seperate out CpG clusters for correlations
Cluster.Ranges <- clusterSitesToGR(clusters)  ##  52707 clusters 12:12 / 467671 clusters 9:12
width(Cluster.Ranges)
mean(width(Cluster.Ranges))   ## 6bp min / 53bp mean / 370bp max  (9:12)

Meth.Clusters <- subsetByOverlaps(cov.lim, granges(clusters))
  head(totalReads(Meth.Clusters))
  head(methReads(Meth.Clusters))
  covStatistics(Meth.Clusters)        ## 159704 CpGs into 22749 clusters (12/12)
                                      ## 250000 - 465000 CpGs into 52707 clusters (9/12) - median 13x

meth.cluster.rel <- rawToRel(Meth.Clusters)
cluster.methylation<-methLevel(meth.cluster.rel)
    colMeans(cluster.methylation,na.rm=T)    #average methylation across samples
    barplot(colMeans(cluster.methylation,na.rm=T),col=c("azure3","azure3","azure3","azure3","azure3","azure3","brown3","brown3","brown3","brown3","brown3","brown3"))     
           
# smooth methylation
smooth <- predictMeth(object = clusters) # h = bandwidth = 80bp
smooth
head(methLevel(smooth))

### Test group effect
ctrl <- smooth[, colData(smooth)$group == "ctrl"]
vitC <- smooth[, colData(smooth)$group == "-vitc"]
mean.ctrl <- rowMeans(methLevel(ctrl))
mean.vitC <- rowMeans(methLevel(vitC))
?plot
plot(mean.ctrl,
       mean.vitC,
       col = "grey",
       xlab = "Methylation in ctrl",
       ylab = "Methylation in -vitc")
      abline(a=0,b=1,col="black",cex=5)
## another way to plot     
      plot(rowMeans(methLevel(smooth[,1:6])),
           rowMeans(methLevel(smooth[,7:12])),
           col = "grey",
           xlab = "Methylation in ctrl",xlim=c(0,1),
           ylab = "Methylation in -vitc",ylim=c(0,1),
           main = "Average methylation - All smoothed clusters")
      abline(a=0,b=1,col="black",cex=5)
## add density 
      dev.off()
      ###load color library
      library(RColorBrewer)
      ###make heat map colors
      Lab.palette <- colorRampPalette(blues9, space = "Lab")
      smoothScatter(mean.ctrl,mean.vitC, cex=1, xlim=c(0,1), ylim=c(0,1),
                  colramp = Lab.palette, nbin=800,bandwidth = .05,nrpoints=Inf)
      
                                  main="RRBS Average Methylation", xlab="Horizontal Location", ylab="Vertical Location")
                  
#-----------------------------------------------------------------------------------------------------------------------------#

# Identify differentially methylated regions (DMRs) 

DMR <- betaRegression(formula = ~group,
                      link = "probit",
                      object = smooth,
                      type = "BR",
                      mc.cores = 4)
head(DMR)

# Parse out significant DMRs 
SigDMR.05 <- findDMRs(DMR,      ### 2158 ranges (DMRs) are found with p < 0.05
                      alpha = 0.05,
                      max.dist = 20,
                      diff.dir = T)

SigDMR.05
quantile(width(SigDMR.05))
mean(width(SigDMR.05))

#-----------------------------------------------------------------------------------------------------------------------------#

## Split into gain and loss of methylation with -VitC 

SigDMR.05
hist(elementMetadata(SigDMR.05)[,4],breaks=100,
     xlim=c(-.4,.4),
     xlab="methylation difference (-VitC - Ctrl)",ylab="number of CpG sites",
     main="Significant methylation sites (p<0.05) E13.5fPGCs")
box()

(moreInVC <- SigDMR.05[elementMetadata(SigDMR.05)[,4]>0])   #1202
(lessInVC <- SigDMR.05[elementMetadata(SigDMR.05)[,4]<0])   #956
	# export(moreInVC, "DMR.05.moreVC.bed")
	# export(lessInVC, "DMR.05.lessVC.bed")
	# export(subsetByOverlaps(covered_regions, SigDMR.05), "DMRs.all.9of12.bed")
	# export(subsetByOverlaps(covered_regions, moreInVC), "DMR_moreVC.Overlap.bed")
	# export(subsetByOverlaps(covered_regions, lessInVC), "DMR_lessVC.Overlap.bed")

elementMetadata(SigDMR.05)
(DMR_5percent_granges <- SigDMR.05[abs(elementMetadata(SigDMR.05)[,4]) > 0.05])
(moreInVC_5percent <- SigDMR.05[elementMetadata(SigDMR.05)[,4]>0.05])     #285 ranges with >5% more methylation in VC removal
(lessInVC_5percent <- SigDMR.05[elementMetadata(SigDMR.05)[,4]<(-0.05)])  #175 ranges with >5% less methylation in VC removal
     hist(elementMetadata(moreInVC_5percent)[,4],breaks=30,
          xlim=c(-.4,.4),ylim=c(0,60),
          xlab="methylation difference (-VitC - Ctrl)",ylab="number of CpG sites",
          main="Significant methylation sites (p<0.05) E13.5fPGCs",sub="5% minimum")
     hist(elementMetadata(lessInVC_5percent)[,4], breaks=15, add=T)
     box()

covered_regions <- clusterSitesToGR(clusters)
covered_regions
quantile(width(covered_regions))
hyper5 <- subsetByOverlaps(covered_regions,moreInVC_5percent) 
quantile(width(hyper5))
hypo5 <- subsetByOverlaps(covered_regions,lessInVC_5percent)
quantile(width(hypo5))

#-----------------------------------------------------------------------------------------------------------------------------#

# FOR MOTIF ANALYSIS NEED ALL GRANGES TO BE 500BP - use flank to increase each bed file 
background_500 <- flank(covered_regions, 500)
quantile(width(background_500))
hyper_500 <- flank(hyper5,500)
hypo_500 <- flank(hypo5,500)
# Export bed files for motif analysis
# export(background_500, "background_9.12_500bp.bed")
# export(hyper_500, "hyperDMR_5percent_500bp.bed")
# export(hypo_500, "hypoDMR_5percent_500bp.bed")
# 
# export(covered_regions,"background_9.12_BroadRRBS.bed")
# export(subsetByOverlaps(covered_regions,moreInVC_5percent), "hyperDMR_5percent_BroadRRBS.bed")
# export(subsetByOverlaps(covered_regions,lessInVC_5percent), "hypoDMR_5percent_BroadRRBS.bed")

#-----------------------------------------------------------------------------------------------------------------------------#

## PLOT HYPER-METH
(percent5_hyper <- subsetByOverlaps(smooth, moreInVC_5percent))
covStatistics(percent5_hyper)        
meth.5percent.hyper <- (methLevel(percent5_hyper))  
barplot(colMeans(meth.5percent.hyper,na.rm=T),col=c("azure3","azure3","azure3","azure3","azure3","azure3","brown3","brown3","brown3","brown3","brown3","brown3"))
#CpG average methylation (ctrl vs -vitc)
head(rowMeans(meth.5percent.hyper[,1:6]))
head(rowMeans(meth.5percent.hyper[,7:12]))
plot(rowMeans(meth.5percent.hyper[,1:6]),
     rowMeans(meth.5percent.hyper[,7:12]),
     col = "red",
     xlim=c(0,1),ylim=c(0,1),
     xlab = "Methylation in ctrl",
     ylab = "Methylation in -vitc",
     main = "Average methylation of CpGs in >5% hypermeth Clusters")
abline(a=0,b=1,col="black",cex=5)

## PLOT HYPO-METH
(percent5_hypo <- subsetByOverlaps(smooth, lessInVC_5percent))
covStatistics(percent5_hypo)        
meth.5percent.hypo <- (methLevel(percent5_hypo))  
barplot(colMeans(meth.5percent.hypo,na.rm=T),col=c("azure3","azure3","azure3","azure3","azure3","azure3","brown3","brown3","brown3","brown3","brown3","brown3"))
#CpG average methylation (ctrl vs -vitc)
plot(rowMeans(meth.5percent.hypo[,1:6]),
     rowMeans(meth.5percent.hypo[,7:12]),
     col = "red",
     xlim=c(0,1),ylim=c(0,1),
     xlab = "Methylation in ctrl",
     ylab = "Methylation in -vitc",
     main = "Average methylation of CpGs in >5% hypermeth Clusters")
abline(a=0,b=1,col="black",cex=5)

#-----------------------------------------------------------------------------------------------------------------------------#

(moreInVC_10percent <- SigDMR.05[elementMetadata(SigDMR.05)[,4]>0.10])     #100 ranges with >10% more methylation in VC removal
(lessInVC_10percent <- SigDMR.05[elementMetadata(SigDMR.05)[,4]<(-0.10)])  #28 ranges with >10% less methylation in VC removal
export(subsetByOverlaps(covered_regions,moreInVC_10percent), "hyperDMR_10percent_BroadRRBS.bed")
export(subsetByOverlaps(covered_regions,lessInVC_10percent), "hypoDMR_10percent_BroadRRBS.bed")

#-----------------------------------------------------------------------------------------------------------------------------#
      
## compare RRBS with RNAseq      
        
# granges of RNAseq differentially expressed genes
up.granges
upgene.meth <- subsetByOverlaps(smooth, up.granges)
  covStatistics(upgene.meth)
  methLevel(upgene.meth)
  boxplot(methLevel(upgene.meth),col=colorBar,main="UpRegulated Gene methylation")
  dev.off()
  heatmap.2(methLevel(upgene.meth),scale="none", col=green2red, cexRow=.1,na.rm=T,
            trace="none",density.info="none",ColSideColors=colorBar)
  plot(rowMeans(methLevel(upgene.meth[,1:6])),
       rowMeans(methLevel(upgene.meth[,7:12])),
       col = "blue",
       xlim=c(0,1),ylim=c(0,1),
       xlab = "Methylation in ctrl",
       ylab = "Methylation in -vitc",
       main = "Average methylation of CpGs across Upregulated Genes")
  abline(a=0,b=1,col="black",cex=5)
  barplot(colMeans(methLevel(upgene.meth),na.rm=T),col=c("azure3","azure3","azure3","azure3","azure3","azure3","brown3","brown3","brown3","brown3","brown3","brown3"))
  box()   

down.granges
downgene.meth <- subsetByOverlaps(smooth, down.granges)
  covStatistics(downgene.meth)
  methLevel(downgene.meth)
  boxplot(methLevel(downgene.meth),col=colorBar,main="downRegulated Gene methylation")
  dev.off()
  heatmap.2(methLevel(downgene.meth),scale="none", col=green2red, cexRow=.1,na.rm=T,
            trace="none",density.info="none",ColSideColors=colorBar,
            dendrogram='column', Rowv=FALSE, Colv=T)
  plot(rowMeans(methLevel(downgene.meth[,1:6])),
       rowMeans(methLevel(downgene.meth[,7:12])),
       col = "blue",
       xlim=c(0,1),ylim=c(0,1),
       xlab = "Methylation in ctrl",
       ylab = "Methylation in -vitc",
       main = "Average methylation of CpGs across Downregulated Genes")
  abline(a=0,b=1,col="black",cex=5)        
  barplot(colMeans(methLevel(downgene.meth),na.rm=T),col=c("azure3","azure3","azure3","azure3","azure3","azure3","brown3","brown3","brown3","brown3","brown3","brown3"))

  