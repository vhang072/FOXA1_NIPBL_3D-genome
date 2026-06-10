determine.overlap.peaks.bed <- function(rep1,rep2,min.overlap.frac){
  library(IRanges)
  chroms <- intersect(unique(rep1[,1]),unique(rep2[,1]))
  total.info1 <- NULL
  total.info2 <- NULL
  for (sp in chroms){
    rep1.chrom <- rep1[rep1[,1] == sp,]
    rep2.chrom <- rep2[rep2[,1] == sp,]
    rep1.overlap.sign <- rep(NA,nrow(rep1.chrom))
    rep2.overlap.sign <- rep(NA,nrow(rep2.chrom))
    for (i in 1:nrow(rep1.chrom)) {
      peak.check.start <- rep1.chrom[i,2]
      peak.check.end <- rep1.chrom[i,3]
      peaks.tomatched.start <- rep2.chrom[,2]
      peaks.tomatched.end <- rep2.chrom[,3]
      f1 <- (peak.check.start <= peaks.tomatched.start & peak.check.end >= peaks.tomatched.start) | 
        (peak.check.start>= peaks.tomatched.start & peak.check.start <= peaks.tomatched.end)
      if (sum(f1)>=1){
        overlap.lengths <- c()
        for (j in 1:sum(f1)){
          t.sort <- sort(c(peak.check.end, peak.check.start,peaks.tomatched.end[f1][j],peaks.tomatched.start[f1][j]))
          overlap.lengths[j] <- (t.sort[3] - t.sort[2] + 1)/min(peak.check.end-peak.check.start+1,peaks.tomatched.end[f1][j]-peaks.tomatched.start[f1][j]+1)
          if (overlap.lengths[j] >= min.overlap.frac){
            rep1.overlap.sign[i] <- "Overlap"
            t.index <- findMatches(T,f1)@to[j]
            rep2.overlap.sign[t.index] <- "Overlap"
          }
        }
      }
    }
    total.info1 <- rbind(total.info1,cbind(rep1.chrom,rep1.overlap.sign))
    total.info2 <- rbind(total.info2,cbind(rep2.chrom,rep2.overlap.sign))
  }
  total.info1 <- as.data.frame(total.info1)
  colnames(total.info1)[ncol(total.info1)] <- "IfOverlap" 
  total.info1$IfOverlap[is.na(total.info1$IfOverlap)==T] <- "NotOverlap"
  total.info2 <- as.data.frame(total.info2)
  colnames(total.info2)[ncol(total.info2)] <- "IfOverlap" 
  total.info2$IfOverlap[is.na(total.info2$IfOverlap)==T] <- "NotOverlap"
  return(list(sample1 = total.info1, sample2 = total.info2))
}


nipbl_peaks <- read.table("./NIPBL.narrowPeak",
                          sep = "\t",header = F)
foxa1_peaks <- read.table("./FOXA1.narrowPeak",
                          sep = "\t",header = F)
foxa1_peaks$V1 <- paste0("chr",foxa1_peaks$V1)
overlaps <- determine.overlap.peaks.bed(nipbl_peaks,foxa1_peaks,0.6)
common <- overlaps$sample1[overlaps$sample1$IfOverlap == "Overlap",]
nipblspe <- overlaps$sample1[overlaps$sample1$IfOverlap == "NotOverlap",]
foxa1spe <- overlaps$sample2[overlaps$sample2$IfOverlap == "NotOverlap",]

write.table(common,"./common.bed",
            sep = "\t",col.names = F, row.names = F, quote = F)
write.table(nipblspe,"./NIPBLspe.bed",
            sep = "\t",col.names = F, row.names = F, quote = F)
write.table(foxa1spe,"./FOXA1spe.bed",
            sep = "\t",col.names = F, row.names = F, quote = F)
