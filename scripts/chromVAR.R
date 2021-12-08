library("chromVAR")
library("motifmatchr")
library("Matrix")
library("SummarizedExperiment")
library("BiocParallel")
library("BSgenome.Hsapiens.UCSC.hg38")
library("tidyverse")

set.seed(2017)
BiocParallel::register(BiocParallel::MulticoreParam(16,progressbar = TRUE))

mount="~/Ambs_ATACseq"

#peaks
# peakfile <- system.file("extdata/test_bed.txt", package = "chromVAR")
peakfile=paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/fixed_width_peaks/regions_of_interest.m1.narrowPeak.annotations.DEGs.Design_3.2.significant.bed",sep="/")
peaks <- getPeaks(peakfile, sort_peaks = TRUE)
peaks

#count
# bamfile <- system.file("extdata/test_RG.bam", package = "chromVAR")
bamfile=paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/qsortedBam/HCC70_1.qsorted.bam",sep="/")
fragment_counts <- getCounts(bamfile, 
                             peaks, 
                             paired =  TRUE, 
                             by_rg = FALSE, 
                             format = "bam",
                             colData = DataFrame(celltype = c("AA")))
colSums(assay(fragment_counts))                         

bamfile1=paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HCC70_1.tn5nicks.bam",sep="/")
fragment_counts1 <- getCounts(bamfile1, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bam",
                             colData = DataFrame(celltype = c("AA")))
colSums(assay(fragment_counts1))


bamfile = c(paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HS578t_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HS578t_2.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/MDAMB436_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/MDAMB436_2.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/BT549_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/BT549_2.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/MDAMB231_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/MDAMB231_2.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HCC1806_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HCC1806_2.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HCC2157_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HCC2157_2.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/MDAMB157_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/MDAMB157_2.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HCC70_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HCC70_2.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/MDAMB468_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/MDAMB468_2.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HCC1500_1.tn5nicks.bam",sep="/"),
            paste(mount,"analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/tn5nicks/HCC1500_2.tn5nicks.bam",sep="/"))
bamfile[2]

# fragment_counts are actually tn5 nick counts
fragment_counts <- getCounts(bamfile, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bam", 
                             colData = DataFrame(celltype = c(rep("AA",8),rep("EA",12))))

colSums(assay(fragment_counts))
fragment_counts

# verify there are no all zero rows
table(rowSums(assay(fragment_counts))==0)

save.image(file="~/Ambs_ATACseq/analysis/project1/CCBR_ATACseq_102621/results/dedupBam/RData")

rowData(fragment_counts)
fragment_counts <- addGCBias(fragment_counts,
                             genome=BSgenome.Hsapiens.UCSC.hg38)
rowData(fragment_counts)
hist(rowData(fragment_counts)$bias)

# create background peaks with the same GC bias
bground <- getBackgroundPeaks(object = fragment_counts)

# colnames(rowData(example_counts))
# #find indices of samples to keep
# counts_filtered <- filterSamples(example_counts, min_depth = 1500, 
#                                  min_in_peaks = 0.15, shiny = FALSE)
# counts_filtered <- filterPeaks(counts_filtered, non_overlapping = TRUE)
motifs <- getJasparMotifs()

# find motifs in ROI
motif_ix <- matchMotifs(motifs, fragment_counts, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)

# get deviations
dev <- computeDeviations(object = fragment_counts, annotations = motif_ix, background_peaks=bground)

# get variations
variability <- computeVariability(dev)

plotVariability(variability, use_plotly = FALSE ) 

vdf=as.data.frame(variability)
devzdf=as.data.frame(assays(dev)$z)
rownames_to_column(vdf,var="motif")-> vdf
vdf=vdf[order(vdf$p_value_adj,-vdf$variability),]
rownames_to_column(devzdf,var="motif")-> devzdf
df=merge(vdf,devzdf,by="motif")
colnames(df)
colnames(df)=gsub(pattern = ".tn5nicks.bam",replacement = "",colnames(df))
colnames(df)

write.table(df,file="chromVAR.results.tsv",row.names=FALSE,col.names = TRUE,quote = FALSE,sep="\t")

topmotifs=head(vdf$motif,100)
topdevzdf=devzdf[(devzdf$motif %in% topmotifs),]
rownames(topdevzdf)=NULL
column_to_rownames(as.data.frame(topdevzdf),var="motif")-> topdevzdf
colnames(topdevzdf)=gsub(pattern = ".tn5nicks.bam",replacement = "",colnames(topdevzdf))
heatmap(as.matrix(topdevzdf))


