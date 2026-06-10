#integrating all samples enrichment
#loading files output from A1 script
data1 <- read.table("./SacCer3.class.txt",
                    sep = "\t",header = T)
data2 <- read.table("./dm6.class.txt",
                    sep = "\t",header = T)
data3 <- read.table("./mm10_liver.class.txt",
                    sep = "\t",header = T)
data4 <- read.table("./mm10_lung.class.txt",
                    sep = "\t",header = T)
data5 <- read.table("./mm10_stomach.class.txt",
                    sep = "\t",header = T)
data6 <- read.table("./hg19_22Rv1.class.txt",
                    sep = "\t",header = T)
data7 <- read.table("./hg19_HCT116.class.txt",
                    sep = "\t",header = T)
data8 <- read.table("./hg19_LoVo.class.txt",
                    sep = "\t",header = T)

sample_name <- c(rep("SacCer",nrow(data1)),
                rep("DM",nrow(data2)),
                rep("Mus_Liver",nrow(data3)),
                rep("Mus_Lung",nrow(data4)),
                rep("Mus_Stomach",nrow(data5)),
                rep("Homo_22Rv1",nrow(data6)),
                rep("Homo_HCT116",nrow(data7)),
                rep("Homo_LoVo",nrow(data8)))
merged <- rbind(data1,data2,data3,data4,data5,data6,data7,data8)
merged$sample_name <- sample_name
merged$signf <- -log10(merged$pval + 1e-100)

merged_f <- merged[merged$fg_frac > 0.15 | merged$signf > 20 ,]

merged_used <- merged[merged$classID %in% unique(merged_f$classID),]

used_motif <- unique(merged_f$classID)
total_occupy <- NULL
for (idx in 1:length(used_motif)){
  total_occupy <- c(total_occupy,sum(merged_used$fg_frac[merged_used$classID %in% used_motif[idx]]))
}
sort_motif <- used_motif[order(total_occupy,decreasing = F)]
merged_used$sample_name <- factor(merged_used$sample_name,levels = c("SacCer","DM","Mus_Liver","Mus_Lung","Mus_Stomach",
                                                                     "Homo_22Rv1","Homo_HCT116","Homo_LoVo"))
merged_used$classID <- factor(merged_used$classID ,levels = sort_motif)

library(ggplot2)
g1 <- ggplot(merged_used,aes(x=sample_name,y=classID,size = fg_frac,fill = signf)) + geom_point(shape = 21) +
  scale_fill_gradientn(colors = rev(rainbow(10,start = 0,end = 2/3)))+ theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
pdf("./enriched_motif_classes_improved.pdf",height = 10,width = 9)
print(g1)
dev.off()
