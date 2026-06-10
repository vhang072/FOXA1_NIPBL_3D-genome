
ctcf <- read.csv("I:/Projects/P2.FOXA1-3D-genome/D0.HiC-analyses_#FOXA1/ChIPseq/ChIP_CTCF_dTAG_0h_SPMR_blacklist.narrowPeak",sep = "\t", header = F)
ctcf[,1] <- paste0("chr",ctcf[,1])
ctcf.sites <- ctcf[ctcf[,5]>50,]

nipbl_spe <- read.table("I:/Projects/P2.FOXA1-3D-genome/D0.HiC-analyses_#FOXA1/ChIPseq/NIPBL-FOXA1.peaks/only.NIPBL.peak.bed",sep = "\t",header = F)
nipbl_common <- read.table("I:/Projects/P2.FOXA1-3D-genome/D0.HiC-analyses_#FOXA1/ChIPseq/NIPBL-FOXA1.peaks/both.NIPBL-FOXA1.peak.bed",sep = "\t",header = F)

final.nipbl <- rbind(nipbl_common,nipbl_spe)
final.nipbl$type <- c(rep("Common",nrow(nipbl_common)),
				rep("NIPBL_spe",nrow(nipbl_spe)))


#novel nipbl appear where
final.nipbl.add.min.dis <- data.frame()
idx <- 0
chroms <- unique(final.nipbl$V1)
for (chrom in chroms){
  print(chrom)
  final.nipbl.chrom <- final.nipbl[final.nipbl[,1] ==chrom,]
  ctcf.sites.chrom <- ctcf.sites[ctcf.sites[,1]==chrom,]
  mid.pos <- as.numeric(apply(as.matrix(ctcf.sites.chrom[,2:3]),1,mean))
  for (i in 1:nrow(final.nipbl.chrom)) {
    if (ceiling(i/100)*100 == i){
      print(i)
    }
    tmp.sort <- sort(abs( mid.pos - mean(as.numeric(final.nipbl.chrom[i,2:3]))),decreasing = F) 
    final.nipbl.add.min.dis[idx,1:5] <- c(final.nipbl.chrom[i,],tmp.sort[1])
    idx <- idx + 1
  }
}
colnames(final.nipbl.add.min.dis)[5] <- "min.dis.to.ctcf.sites"
final.nipbl.add.min.dis[,5] <- as.numeric(final.nipbl.add.min.dis[,5])
final.nipbl.add.min.dis$type <- as.factor(final.nipbl.add.min.dis$type)
library(ggplot2)
pdf("./nipbl2ctcf.pdf",height = 3.5,width = 5)
g <- ggplot(final.nipbl.add.min.dis, aes(x=log10(min.dis.to.ctcf.sites+1),  color = type)) +  xlim(0,7) + geom_density(alpha =0) +
  theme_bw()
print(g)
dev.off()


