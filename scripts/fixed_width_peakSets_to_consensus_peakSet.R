#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-i", "--inputNarrowPeaks", 
                    type="character", 
                    help="fixed width narrowPeak input files from narrowPeak2uniformwidthpeakset.R",
                    required=TRUE)
parser$add_argument("-o", "--outputNarrowPeak", 
                    type="character", 
                    help="narrowPeak output file",
                    required=TRUE)
parser$add_argument("-t", "--tmpdir", 
                    type="character",
                    help="tmp dir",
                    required=FALSE,
                    default=NULL)

args <- parser$parse_args()

norm_pvalue_threshold=5
min_replicates=2

suppressPackageStartupMessages(library("tidyverse"))

debug=0

narrowPeaks=unlist(strsplit(args$inputNarrowPeak,","))
out_narrowPeak=args$outputNarrowPeak

if (debug==1){
narrowPeaks=unlist(strsplit("/Volumes/Ambs_ATACseq/analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/uniform_width_peaks/HCC2157_1.genrich.samplePeakSet.narrowPeak,/Volumes/Ambs_ATACseq/analysis/project1/CCBR_ATACseq_102621/results/peaks/genrich/uniform_width_peaks/HCC2157_2.genrich.samplePeakSet.narrowPeak",","))
out_narrowPeak="HCC2157.genrich.consensus.narrowPeak"
setwd(dirname(narrowPeak))
} else {
  if(is.null(args$tmpdir)){
    tmpdir=setwd(dirname(out_narrowPeak))
  } else {
    tmpdir=args$tmpdir
  }  
  setwd(tmpdir)
}
bn=basename(out_narrowPeak)

readnarrowpeak<-function(narrowPeak){
  # cat(narrowPeak)
  bn=basename(narrowPeak)
  x=read.csv(narrowPeak,
             header=FALSE,
             sep = "\t",
             check.names = FALSE,
             strip.white = TRUE
  )

  colnames(x)=c("chrom",
                "start",
                "end",
                "peakname",
                "score",
                "strand",
                "signalValue",
                "normalizedpvalue",
                "qvalue",
                "summitdist")
  x$peakname=paste(bn,x$peakname,sep="_")
  return(x)
}

xdf <- readnarrowpeak(narrowPeaks[1])

for (i in 2:length(narrowPeaks)) {
  xdf2=readnarrowpeak(narrowPeaks[i])
  xdf=rbind(xdf,xdf2)
}

tmp1bed=paste(bn,"tmp1.bed",sep=".")
tmp=paste(bn,"tmp",sep=".")
mutate_at(xdf,vars("start","end","summitdist"),as.integer) -> xdf
arrange(xdf,desc(normalizedpvalue)) -> xdf
xdf %>% mutate(rank=row_number()) -> xdf

write.table(xdf,file=tmp1bed,
            sep="\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
# q()
cmd=paste("bedtools intersect -wa -wb -a",tmp1bed,"-b",tmp1bed,">",tmp)
print(cmd)
system(cmd)
x2=read.csv(tmp,
           header=FALSE,
           sep = "\t",
           check.names = FALSE,
           strip.white = TRUE
)
colnames(x2)=c("chrom",
              "start",
              "end",
              "peakname",
              "score",
              "strand",
              "signalValue",
              "normalizedpvalue",
              "qvalue",
              "summitdist",
              "rank",
              "chrom2",
              "start2",
              "end2",
              "peakname2",
              "score2",
              "strand2",
              "signalValue2",
              "normalizedpvalue2",
              "qvalue2",
              "summitdist2",
              "rankb")


keep1_ranks=x2[x2$rank==x2$rankb,]$rank
remove_ranks=x2[x2$rank<x2$rankb,]$rankb
keep_ranks=setdiff(keep1_ranks,remove_ranks)

xdf2=xdf[xdf$rank %in% keep_ranks,]

xdf2=xdf2[,c("chrom",
                  "start",
                  "end",
                  "peakname",
                  "score",
                  "strand",
                  "signalValue",
                  "normalizedpvalue",
                  "qvalue",
                  "summitdist")]
tmp2bed=paste(bn,"tmp2.bed",sep=".")
mutate_at(xdf2,vars("start","end","summitdist"),as.integer) -> xdf2
# write merged peaks
write.table(xdf2,file=tmp2bed,
            sep="\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


peaksoverlapcount=data.frame(peakname=xdf2$peakname,count=0)
for (i in 1:length(narrowPeaks)) {
  cmd=paste("bedtools intersect -wa -wb -a",tmp2bed,"-b",narrowPeaks[i],">",tmp)
  print(cmd)
  system(cmd)
  x2=read.csv(tmp,
              header=FALSE,
              sep = "\t",
              check.names = FALSE,
              strip.white = TRUE
  )
  colnames(x2)=c("chrom",
                 "start",
                 "end",
                 "peakname",
                 "score",
                 "strand",
                 "signalValue",
                 "normalizedpvalue",
                 "qvalue",
                 "summitdist",
                 "chrom2",
                 "start2",
                 "end2",
                 "peakname2",
                 "score2",
                 "strand2",
                 "signalValue2",
                 "normalizedpvalue2",
                 "qvalue2",
                 "summitdist2")
  x2[x2$normalizedpvalue2>norm_pvalue_threshold,]->x2
  update_peaks=peaksoverlapcount$peakname %in% x2$peakname
  peaksoverlapcount[update_peaks, "count"] <- 1 + peaksoverlapcount[update_peaks, "count"]
  
  system(paste("rm -f",tmp))
}
keep_peaks=peaksoverlapcount[peaksoverlapcount$count >= min_replicates,]$peakname
xdf3=xdf2[xdf2$peakname %in% keep_peaks,]
write.table(xdf3,file=out_narrowPeak,
            sep="\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
system(paste("rm -f",tmp1bed,tmp2bed))
