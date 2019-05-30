if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPseeker")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("affy")
BiocManager::install("estrogen")
BiocManager::install("vsn")
BiocManager::install("genefilter")

library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
install.packages("ggplot2")
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)
library(ChIPseeker)

BiocManager::install("DiffBind")
library(DiffBind)

#loading sample data
setwd(system.file("extra", package="DiffBind"))
list.files()
read.csv("tamoxifen.csv")

#calculate differential binding analysis 
tamoxifen <- dba(sampleSheet="tamoxifen.csv")

#plot correlation heatmap based on ER occupancy
plot(tamoxifen)

#Counting reads
#tamoxifen <- dba.count(tamoxifen, summits=250)
tamoxifen2 <- data("tamoxifen_counts") #gives us the number of reads per peak as a value of "strength" or certainty of bidning

## COMPARE DIFFERENTIAL BINDING OF ER BASED ON SENSIVIE VS. RESISTANT CELL LINES
#Establishing the condition of which to sort (column 4 of the tamoxifen.csv file)
tamoxifen2 <- dba.contrast(tamoxifen, categories=DBA_CONDITION)
tamoxifen2 <- dba.analyze(tamoxifen2) #Runs DESEQ2 analysis


#Visualize with heatmap
dev.off()
plot(tamoxifen2, contrast = 1)
#correlation heatmap based only on DIFFERNTIALLY BOUND regions

#Retreiving differntially bound sites
tamoxifen.DB <- dba.report(tamoxifen2)
tamoxifen.DB  #Returns a GRanges object

#The value columns show the mean read concentration over all the samples (log2 normalized ChIP read counts with control read counts subtracted) 
#and the mean concentration over the first (Resistant) group and second (Responsive) group. 
#Also gives FOLD change which is indicated by positive values showing increase affinity and negative values as decreased
#P-value and corrected FDR

#VISUALIZING DIFFERENTIAL BINDING
#PCA Analysis
data(tamoxifen_analysis) 
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_CONDITION)
dba.plotPCA(tamoxifen, contrast=1,label=DBA_TISSUE)

#Volcano plots
dba.plotVolcano(tamoxifen2)

#Box plots
sum(tamoxifen.DB$Fold<0)
sum(tamoxifen.DB$Fold>0)
pvals <- dba.plotBox(tamoxifen)

#Heatmaps
corvals <- dba.plotHeatmap(tamoxifen)
corvals <- dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE)


#vizualization of features of ChIP binding profile
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
#Read the ER ChIP profiles for Tamoxifen sensitive and resistant MCF7 samples
list.files("peaks")
MCF7 <-readPeakFile("peaks/MCF7_ER_1.bed.gz", as = "GRanges")
Resistant <-readPeakFile("peaks/TAMR_ER_1.bed.gz", as ="GRanges")

#Look at the ChIP binding peaks throughout chromosome 18
dev.off() #fixes a weird issue with ggplot
covplot(MCF7, weightCol="V5")
#covplot(MCF7, weightCol="V5", chrs=c("chr18"), xlim=c(4.5e7, 5e7))
covplot(Resistant, weightCol="V5")
#covplot(Resistant, weightCol="V5", chrs=c("chr18"), xlim=c(4.5e7, 5e7))


###PROMOTER BINDING PROFILE
#get binding within promoter regions: SENSITIVE
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_MCF7 <- getTagMatrix(MCF7, windows=promoter)
tagHeatmap(tagMatrix_MCF7, xlim=c(-3000, 3000), color="red")

#get binding within promoter regions: RESISTANT
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_resistant <- getTagMatrix(Resistant, windows=promoter)
tagHeatmap(tagMatrix_resistant, xlim=c(-3000, 3000), color="blue")

#binding within promoter regions: RESISTANT, TRIPLE NEGATIVE BREAST CANCER
BT474 <-readPeakFile("peaks/BT474_ER_1.bed.gz", as ="GRanges")
tagMatrix_BT474 <- getTagMatrix(BT474, windows=promoter)
tagHeatmap(tagMatrix_BT474, xlim=c(-3000, 3000), color="green")

###AVERAGE READ COUNT FREQUENCY FOR PROMOTER
#MCF7
plotAvgProf(tagMatrix_MCF7, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
#RESISTANT
plotAvgProf(tagMatrix_resistant, xlim=c(-3000, 3000),
        xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
#BT474
plotAvgProf(tagMatrix_BT474, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

##ANNOTATING BINDING REGIONS 
##Identification of the functional annotation of the region bound by the protein of interest
MCF7_peakAnno <- annotatePeak(MCF7, tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Hs.eg.db")
dev.off() #another weird ggplot error
plotAnnoPie(MCF7_peakAnno)

#other plots(probably less helpful)
plotAnnoBar(MCF7_peakAnno)
vennpie(MCF7_peakAnno)
upsetplot(MCF7_peakAnno)

plotDistToTSS(MCF7_peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

#What genes are bound by the TF in the promoter region?

library(ReactomePA)
pathway1 <- enrichPathway(as.data.frame(MCF7_peakAnno)$geneId)
head(pathway1, 2)
dotplot(pathway1)

#gene <- seq2gene(MCF7, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
#pathway2 <- enrichPathway(gene)
#head(pathway2, 2)


##WHAT ABOUT THE TAMOXIFEN RESISTANT CELLS
#ERres_peakAnno <- annotatePeak(Resistant, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
#plotAnnoPie(ERres_peakAnno)
#pathway_res <- enrichPathway(as.data.frame(ERres_peakAnno)$geneId)
#head(pathway_res,2)
#dotplot(pathway_res)


#Integrated analysis with RNA data
#New example
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

files <- getSampleFiles()
print(files)

peak <- readPeakFile(files[[4]])
peak

# to speed up the compilation of this vigenette, we load a precaculated tagMatrixList
data("tagMatrixList")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))


tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

p <- GRanges(seqnames=c("chr1", "chr3"),
             ranges=IRanges(start=c(1, 100), end=c(50, 130)))
shuffle(p, TxDb=txdb)

#Finding overlapping peaks
enrichPeakOverlap(queryPeak     = files[[5]],
                  targetPeak    = unlist(files[1:4]),
                  TxDb          = txdb,
                  pAdjustMethod = "BH",
                  nShuffle      = 50,
                  chainFile     = NULL,
                  verbose       = FALSE)



