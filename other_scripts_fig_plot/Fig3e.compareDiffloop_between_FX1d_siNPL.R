diff_loop1 <- read.table("diff_loops_0hVS24h_improved.txt",
                         sep = "\t",header = T)
diff_loop2 <- read.table("diff_loops_siCtrlVSsiNIPBL_improved.txt",
                         sep = "\t",header = T)

rownames(diff_loop1) <- paste(diff_loop1[,1],diff_loop1[,2],diff_loop1[,5],sep = "_")
rownames(diff_loop2) <- paste(diff_loop2[,1],diff_loop2[,2],diff_loop2[,5],sep = "_")

loopused <- intersect(rownames(diff_loop1),rownames(diff_loop2))
diff1 <- diff_loop1[loopused,]
diff2 <- diff_loop2[loopused,]

mergeinfo <- cbind(diff1[,20],diff2[,20])

tb <- table(mergeinfo[,1],mergeinfo[,2])
data <- matrix(as.numeric(tb),ncol = 3,byrow = F)

rownames(data) <- c("De.FOXA1","In.FOXA1","Unch.FOXA1")
colnames(data) <- c("De.NIPBL","In.NIPBL","Unch.NIPBL")

library(pheatmap)
pdf("./CompareDiffLoops_FX1d_siPBL.pdf",height = 4,width = 5.2)
pheatmap(data,cluster_rows = F, cluster_cols = F,display_numbers = T,number_format = "%d",fontsize_number =15)
dev.off()

mergeddata <- cbind(diff1,diff2[,20])
mergeddata$labelmerge <- paste0(mergeddata[,20],mergeddata[,21])

f1 <- mergeddata[mergeddata$labelmerge %in% c("DecreasedDecreased","DecreasedUnchanged","UnchangedUnchanged","UnchangedDecreased","IncreasedUnchanged",
                                              "UnchangedIncreased"),]

write.table(f1,"./diff_loops_FOXA1d_siNIPBL_combined.txt",
            sep = "\t",row.names = F,quote = F)
