#NIPBL peaks centered at peak summit
data <- read.table("./Trimed200bp_hg19_LoVo.bed",
                   sep = "\t",header = F)
data$mid_pos <- apply(data[,2:3],1,mean)

#TAD boundary file
bound <- read.table("./h0_1_20000_iced.boundary.txt",
                    sep = "\t",header =T)
bound$mid <- floor(apply(bound[,2:3],1,mean))
chroms <- unique(bound[,1])

dis_info <- NULL
for (chrom in chroms){
  data_chrom <- data[data[,1] == chrom,]
  bound_chrom <- bound[bound[,1] == chrom,]
  if (nrow(data_chrom)>=1){
    for (idx in 1:nrow(data_chrom)){
      min_dis <- min(abs(data_chrom[idx,"mid_pos"] -bound_chrom[,"mid"]))
      dis_info <- rbind(dis_info,c(data_chrom[idx,],min_dis))
    }
  }
}


wanted <- dis_info[,c(4,12)]
wanted <- as.data.frame(wanted)
colnames(wanted) <- c("id","mindis2bound")
wanted[,1] <- as.character(wanted[,1])
wanted[,2] <- as.numeric(wanted[,2])

#distribution of distance to TAD boudaries from NIPBL peaks
write.table(wanted,"./hg19_LoVo.txt",
            sep = "\t",col.names = F, row.names = F, quote = F)

library(ggplot2)

ggplot(wanted,aes(x = log10(mindis2bound+100))) + geom_density() + theme_bw()
