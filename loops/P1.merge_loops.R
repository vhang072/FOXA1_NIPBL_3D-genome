##merge loops 
loop1 <- read.table("downsample_FX_dTAG_0h_10kb.bed",sep = "\t",header = F)
loop2 <- read.table("downsample_FX_dTAG_24h_10kb.bed",sep = "\t",header = F)

loop1$sign.overlap <- 0
loop2$sign.overlap <- 0

overlapped.loops <- NULL
dis.threshold <- 20000
for (i in 1:nrow(loop1)) {
  anchor1.chrom <- loop1[i,1]
  anchor1.start <- loop1[i,2]
  anchor2.chrom <- loop1[i,4]
  anchor2.start <- loop1[i,5]
  
  f1 <- anchor1.chrom == loop2[,1] & anchor2.chrom == loop2[,4] & anchor1.start >= (loop2[,2]-dis.threshold) &
    anchor1.start <= (loop2[,2]+dis.threshold) & anchor2.start >= (loop2[,5]-dis.threshold) & anchor2.start <= (loop2[,5]+dis.threshold)
  
  if (sum(f1)==1) {
    dis.anchor1 <- loop2[f1,2]-anchor1.start
    dis.anchor2 <- loop2[f1,5]-anchor2.start
    overlapped.tmp.loop <- cbind(loop1[i,],dis.anchor1,dis.anchor2)
    overlapped.loops <- rbind(overlapped.loops,overlapped.tmp.loop)
    loop1$sign.overlap[i] <- 1
    loop2$sign.overlap[f1] <- 1
  } else if (sum(f1)>1) {
    print(loop1[i,])
    print(loop2[f1,])
    loop1$sign.overlap[i] <- sum(f1)
    loop2$sign.overlap[f1] <- sum(f1)
  }
}

overlapped.loops2 <- overlapped.loops
for (i in 1:nrow(overlapped.loops2)) {
  if (abs(overlapped.loops2[i,8])==20000) {
    overlapped.loops2[i,2] <- overlapped.loops2[i,2] + overlapped.loops2[i,8]/2
    overlapped.loops2[i,3] <- overlapped.loops2[i,3] + overlapped.loops2[i,8]/2
  }
  if (abs(overlapped.loops2[i,9])==20000) {
    overlapped.loops2[i,5] <- overlapped.loops2[i,5] + overlapped.loops2[i,9]/2
    overlapped.loops2[i,6] <- overlapped.loops2[i,6] + overlapped.loops2[i,9]/2
  }
}
overlapped.loops2[,7] <- 1
merged.loops <- rbind(loop1[loop1[,7]==0,],loop2[loop2[,7]==0,],overlapped.loops2[,1:7])

write.table(merged.loops,file="merged.loops.txt",sep = '\t',row.names = F,col.names = T,quote = F)
