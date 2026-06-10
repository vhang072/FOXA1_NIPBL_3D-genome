
#output from A3 script where deposit the distance to TAD boundary for NIPBL or NIPBL-like binding sites
peakdis <- read.table("./hg19_22Rv1.txt",
                      sep = "\t",header = F)
motifdetail <- read.table("./hg19_22Rv1_fg_motifs.txt",
                          sep = "\t",header = T)

#motif information from Jaspar
info <- read.table("./Vertebrata.info",
                   sep = "\t",header = F)
for (idx in 1:5){
  info[,idx] <- str_trim(info[,idx],side = "right")
}
info[info$V4 == "Tryptophan cluster factors","V4"] <- info[info$V4 == "Tryptophan cluster factors","V5"]
info$V4[info$V4==""] <- info$V2[info$V4==""]


#chose motif with high occurence
motifs <- table(motifdetail$Motif.Name)
used_motifs <- motifs[motifs/nrow(peakdis) > 0.05]


density_res <- NULL
used_classes <- NULL
mean_dis <- NULL
for (x in 1:length(used_motifs)){
  target <- names(used_motifs[x])
  target_id <- unlist(str_split(target,":::"))[1]
  target_class <- info[info$V1 == target_id,"V4"]
  peaks <- motifdetail[motifdetail$Motif.Name == target,"PositionID"]
  dis <- peakdis[peakdis$V1 %in% peaks,2]
  out <- density(log10(dis + 1e4))
  if (mean(dis)>mean(peakdis[,2])){
    type_line <- "far"
  } else {
    type_line <- "near"
  }
  mean_dis <- rbind(mean_dis,c(mean(dis),target,target_class,type_line))
  density_res <- rbind(density_res,cbind(out$x,
                                         out$y,
                                         rep(target,length(out$y)),
                                         rep(target_class,length(out$y)),
                                         rep(type_line,length(out$y))))
  used_classes <- c(used_classes,target_class)
}

out <- density(log10(peakdis[,2]+1e4))
density_res <- rbind(density_res,cbind(out$x,
                                       out$y,
                                       rep("NIPBL",length(out$y)),
                                       rep("NIPBL",length(out$y)),
                                       rep("NIPBL",length(out$y))))

density_res <- as.data.frame(density_res)
density_res[,1] <- as.numeric(density_res[,1])
density_res[,2] <- as.numeric(density_res[,2])
density_res[,3] <- as.character(density_res[,3])
density_res[,4] <- factor(density_res[,4], levels = c(sort(unique(used_classes)),"NIPBL"))
density_res[,5] <- as.character(density_res[,5])
colnames(density_res) <- c("log10dis","density","motif","class","type")

library(ggplot2)
pdf("./hg19_22Rv1_dis.pdf",height = 6,width = 20)
ggplot(density_res,aes(x = log10dis,y = density,group = motif,colour = class)) + 
  geom_line(aes(linetype = type),linewidth =0.75) + 
  scale_color_manual(values = rainbow(length(unique(density_res$class)),start = 0,end = 3/4)) +
  scale_linetype_manual(values= c("dashed", "dotted","solid"))+xlim(3.5,6.5)+theme_bw()
dev.off()

mean_dis <- as.data.frame(mean_dis)
mean_dis[,1] <- as.numeric(mean_dis[,1])
mean_dis[,2] <- as.character(mean_dis[,2])
mean_dis[,3] <- as.character(mean_dis[,3])
mean_dis[,4] <- as.character(mean_dis[,4])
#colnames(mean_dis) <- c("mean_dis","motifID","class","type")
#colnames(mean_dis) <- c("mean_dis","motifID","class","type")
selected_motifs <- c("MA1583.1:::ZFP57","MA1464.1:::ARNT2","MA0139.1:::CTCF","MA0028.2:::ELK1","MA0506.1:::NRF1","MA0763.1:::ETV3",
                     "MA0032.2:::FOXC1","MA1487.1:::FOXE1","MA0901.2:::HOXB13","MA1113.2:::PBX2","MA0602.1:::Arid5a","MA0785.1:::POU2F1",
                     "NIPBL")
selected_density <- density_res[density_res$motif %in% selected_motifs,]
selected_density$motif <- factor(selected_density$motif, levels = selected_motifs)
#selected_density$motif <- paste0(selected_density$motif,"_",selected_density$class)
pdf("./hg19_22Rv1_dis_tops.pdf",height = 8,width = 14)
ggplot(selected_density,aes(x = log10dis,y = density,group = motif,colour = motif)) + 
  geom_line(aes(linetype = type),linewidth =0.75) + 
  scale_color_manual(values = rainbow(length(unique(selected_density$motif)),start = 0,end = 3/4)) +
  scale_linetype_manual(values= c("dashed", "dotted","solid"))+theme_bw() + xlim(3.5,6.5)
dev.off()
#ggplot(mean_dis,aes(x = mean_dis)) + geom_histogram() + geom_density(aes(x = V2),data = peakdis) + theme_bw()

write.table(mean_dis,"./hg19_22Rv1_meandis.txt",
            sep = "\t",col.names = F, row.names = F, quote = F)