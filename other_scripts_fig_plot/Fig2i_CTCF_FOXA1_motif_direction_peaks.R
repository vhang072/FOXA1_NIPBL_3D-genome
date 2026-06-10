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

motif_anno <- function(peaks, motif){
  peaks_anno <- NULL
  chroms_used <- unique(peaks$V1)
  peaks_sorted <- NULL
  for (chrom in chroms_used){
    peaks_chrom <- peaks[peaks$V1 == chrom,]
    motif_chrom <- motif[motif$sequence_name == chrom,]
    for (idx in 1:nrow(peaks_chrom)){
      f1 <- motif_chrom[(motif_chrom$start >= peaks_chrom[idx,2]) & (motif_chrom$stop <= peaks_chrom[idx,3]),]
      motif_num <- nrow(f1)
      if (motif_num >= 1){
        strand_unqiue <- unique(f1$strand)
        if (length(strand_unqiue) >=2){
          motif_strand <- "ambiguous"
        } else {
          motif_strand <- strand_unqiue
        }
      } else {
        motif_strand <- NA
      }    
      peaks_anno <- rbind(peaks_anno, c(motif_num,motif_strand))
    }
    peaks_sorted <- rbind(peaks_sorted, peaks_chrom)
  }
  colnames(peaks_anno) <- c("motif_num","motif_strand")
  peaksout <- cbind(peaks_sorted,peaks_anno)
  return(peaksout)
}


ctcf_peaks <- read.table("./ChIP_CTCF_dTAG_0h_SPMR_blacklist.narrowPeak",
                         sep = "\t",header = F)
ctcf_peaks$V1 <- paste0("chr",ctcf_peaks$V1)
nipbl_peaks <- read.table("./FXdTAG_0h_Rv1_NIPBL_blacklist.narrowPeak",
                         sep = "\t",header = F)
nipbl_peaks$V1 <- paste0("chr",nipbl_peaks$V1)
foxa1_peaks <- read.table(".22Rv1_FOXA1_calledpeaks_hg19.bed",
                          sep = "\t",header = F)

#1 choice
# ctcf_motif <- read.table("./meme_fimo_CTCF_MA0139.1_1e-4_hg19.tsv",
#                          sep = "\t",header = T)
# foxa1_motif <- read.table("./meme_fimo_FOXA1_MA0148.1_1e-4_hg19.tsv",
#                          sep = "\t",header = T)
#2 choice
ctcf_motif0 <- read.table("./ctcf.hg19.bed",
                           sep = "\t",header = F)
ctcf_motif <- ctcf_motif0[ctcf_motif0$V5 >=9, ]
colnames(ctcf_motif)[c(1:3,ncol(ctcf_motif))] <- c("sequence_name","start","stop","strand")
foxa1_motif0 <- read.table("./foxa1.lncap.hg19.bed",
                          sep = "\t",header = F)
foxa1_motif <- foxa1_motif0[foxa1_motif0$V5 >=9, ]
colnames(foxa1_motif)[c(1:3,ncol(foxa1_motif))] <- c("sequence_name","start","stop","strand")



chroms_used <- c(paste0("chr",seq(1,22)),"chrX")
ctcf <- ctcf_peaks[ctcf_peaks$V1 %in% chroms_used,]
nipbl <- nipbl_peaks[nipbl_peaks$V1 %in% chroms_used,]
foxa1 <- foxa1_peaks[foxa1_peaks$V1 %in% chroms_used,]

ctcf_motif <- ctcf_motif[ctcf_motif$sequence_name %in% chroms_used,]
foxa1_motif <- foxa1_motif[foxa1_motif$sequence_name %in% chroms_used,]

ctcf <- motif_anno(ctcf,ctcf_motif)
foxa1 <- motif_anno(foxa1,foxa1_motif)


ovlap1 <- determine.overlap.peaks.bed(ctcf,nipbl,0.3)
ctcf2 <- ovlap1$sample1
ctcf3 <- ctcf2[(ctcf2$IfOverlap == "Overlap") & (is.na(ctcf2$motif_strand) == F),]
ovlap2 <- determine.overlap.peaks.bed(ctcf3[,1:(ncol(ctcf3)-1)],foxa1,0.3)
ctcf_onlynipbl <- ovlap2$sample1[ovlap2$sample1$IfOverlap == "NotOverlap",]

ovlap3 <- determine.overlap.peaks.bed(foxa1,nipbl,0.3)
foxa1_2 <- ovlap3$sample1
foxa1_3 <- foxa1_2[(foxa1_2$IfOverlap == "Overlap") & (is.na(foxa1_2$motif_strand) == F),]
ovlap4 <- determine.overlap.peaks.bed(foxa1_3[,1:(ncol(foxa1_3)-1)],ctcf,0.3)
foxa1_onlynipbl <- ovlap4$sample1[ovlap4$sample1$IfOverlap == "NotOverlap",]

write.table(ctcf_onlynipbl, "./ctcf_onlynipbl_motifori.bed",
            sep = "\t",row.names = F, col.names = T)
write.table(foxa1_onlynipbl, "./foxa1_onlynipbl_motifori.bed",
            sep = "\t",row.names = F, col.names = T)