diff_loop1 <- read.table("./diff_loops_FOXA1d_siNIPBL_combined.txt",
                         sep = "\t",header = T)
diff_loop1 <- diff_loop1[diff_loop1$labelmerge %in% c("DecreasedDecreased","DecreasedUnchanged","UnchangedDecreased"),]
tad <- read.table("./FOXA1_dtag_TAD_types.txt",
                  sep = "\t",header = T)



#calculate the relative pos betwwen loop anchor and TAD boundaries
relative_pos <- NULL
for (idx in 1:nrow(diff_loop1)){
  loop_start <- mean(as.numeric(diff_loop1[idx,2:3]))
  loop_end <- mean(as.numeric(diff_loop1[idx,5:6]))
  loop_mid <- mean(c(loop_start,loop_end))
  f1 <- (tad$chrom == diff_loop1[idx,1]) & (tad$start <= loop_mid) & (tad$xEnd >= loop_mid)
  if (sum(f1)==1){
    tad_start <- tad[f1,2]
    tad_end <- tad[f1,3]
    relative_pos1 <- (loop_start - tad_start)/(tad_end - tad_start)
    relative_pos2 <- (loop_end - tad_start)/(tad_end - tad_start)
  } else {
    relative_pos1 <- NA
    relative_pos2 <- NA
  }
  relative_pos <- rbind(relative_pos,c(relative_pos1,relative_pos2))
}

annotation <- cbind(diff_loop1[,22],relative_pos)
annotation <- as.data.frame(annotation)
colnames(annotation) <- c("loop.label","anchor1","anchor2")

annotation[,2] <- as.numeric(annotation[,2])
annotation[,3] <- as.numeric(annotation[,3])
library(ggplot2)
#ggplot(annotation,aes(x = anchor1,y = anchor2,fill = loop.label)) + geom_point(shape = 21,size = 1.0) + theme_bw() + xlim(-1,1) + ylim(0,2)

x1 <- annotation[,c(1,2)]
x2 <- annotation[,c(1,3)]


colnames(x1) <- c("loop","anchor")
colnames(x2) <- c("loop","anchor")
annotation2 <- rbind(x1,x2)
pdf("./Loopanchor2tad_combined.pdf",height = 4,width = 7)
ggplot(annotation2,aes(x = anchor,colour = loop)) + geom_density(size = 1.25) + xlim(-0.2,1.2) +theme_bw()
  #scale_color_manual(values = c("#6A5ACD","#FF8000","#C0C0C0"))
dev.off()




