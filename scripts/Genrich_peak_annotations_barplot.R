mount="~/Ambs_ATACseq/"
setwd(paste0(mount,"analysis/project1/CCBR_ATACseq_102621/results/tnbc_only_QC"))
library(ggplot2)
library(paletteer)
library(ggthemes)
df=read.csv("Genrich_Peak_Annotations_mqc.csv",header = TRUE,sep = "\t",check.names = FALSE,comment.char = "#")
head(df)
moltendf=as.data.frame(melt(df))
colnames(moltendf)=c("cellline","annotation","npeaks")
head(moltendf)
p=ggplot(moltendf, aes(fill=annotation, y=npeaks, x=cellline)) + 
  geom_bar(position="stack", stat="identity", width = 0.5, ,colour="black")+
  coord_flip()+
  theme_classic()+
  theme(legend.title = element_blank(),legend.text=element_text(size=6),legend.position = "bottom")
p + scale_fill_hue(h = c(180, 300))
p + scale_fill_paletteer_d(ptol_pal()(12))
p + scale_color_viridis(discrete=TRUE, option="viridis")

p
pdf("barplot.pdf",paper="USr")
print(p)
dev.off()
