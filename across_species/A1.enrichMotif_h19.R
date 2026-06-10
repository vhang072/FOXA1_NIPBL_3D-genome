library(stringr)
#fg_line_num <- 22657 #22Rv1
#fg_line_num <- 15355 #HCT116
fg_line_num <- 29740 #LoVo
bg_line_num <- 40000

#load files output from Homer motif finder
#files like
#PositionID	Offset	Sequence	Motif Name	Strand	MotifScore
#./peak/FXdTAG_0h_Rv1_NIPBL_SPMR_peak_26725	-25	CACGTG	MA0004.1:::Arnt	+	7.255608
#./peak/FXdTAG_0h_Rv1_NIPBL_SPMR_peak_26725	-20	CACGTG	MA0004.1:::Arnt	-	7.255608
#./peak/FXdTAG_0h_Rv1_NIPBL_SPMR_peak_26723	83	CACGTG	MA0004.1:::Arnt	+	7.255608
#./peak/FXdTAG_0h_Rv1_NIPBL_SPMR_peak_26723	88	CACGTG	MA0004.1:::Arnt	-	7.255608

fg <- read.table("./hg19_LoVo_fg_motifs.txt",
                 sep = "\t",header = T)
fg$motif.ID <- matrix(unlist(str_split(fg$Motif.Name,":::")),ncol = 2,byrow = T)[,1]
bg <- read.table("./hg19_bg_motifs.txt",
                 sep = "\t",header = T)
bg$motif.ID <- matrix(unlist(str_split(bg$Motif.Name,":::")),ncol = 2,byrow = T)[,1]

#motif annotation information from Jaspar
info <- read.table("./Vertebrata.info",
                   sep = "\t",header = F)
for (idx in 1:5){
  info[,idx] <- str_trim(info[,idx],side = "right")
}
info[info$V4 == "Tryptophan cluster factors","V4"] <- info[info$V4 == "Tryptophan cluster factors","V5"]
info$V4[info$V4==""] <- info$V2[info$V4==""]
unique_class <- unique(info$V4)

enrichout <- NULL
for (idx in 1:length(unique_class)){
  this_class_ids <- info$V1[info$V4 == unique_class[idx]]
  fg_unique_peaks <- unique(fg[fg$motif.ID %in% this_class_ids,1])
  bg_unique_peaks <- unique(bg[bg$motif.ID %in% this_class_ids,1])
  
  tbl <- rbind(c(length(fg_unique_peaks),fg_line_num-length(fg_unique_peaks)),
               c(length(bg_unique_peaks)+1,bg_line_num-length(bg_unique_peaks)))
  #Chisq for enrichment test
  test <- chisq.test(tbl)
  enrichout <- rbind(enrichout,c(unique_class[idx],length(this_class_ids),
                                 length(fg_unique_peaks),fg_line_num,
                                 length(bg_unique_peaks),bg_line_num,
                                 length(fg_unique_peaks)/fg_line_num,
                                 length(bg_unique_peaks)/bg_line_num,
                                 test$p.value))
}

enrichout2 <- NULL
#each motif
for (idx in 1:nrow(info)){
  fg_unique_peaks <- unique(fg[fg$motif.ID %in% info$V1[idx],1])
  bg_unique_peaks <- unique(bg[bg$motif.ID %in% info$V1[idx],1])
  
  tbl <- rbind(c(length(fg_unique_peaks),fg_line_num-length(fg_unique_peaks)),
               c(length(bg_unique_peaks)+1,bg_line_num-length(bg_unique_peaks)))
  test <- chisq.test(tbl)
  enrichout2 <- rbind(enrichout2,c(info$V1[idx],
                                   length(fg_unique_peaks),fg_line_num,
                                   length(bg_unique_peaks),bg_line_num,
                                   length(fg_unique_peaks)/fg_line_num,
                                   length(bg_unique_peaks)/bg_line_num,
                                   test$p.value,info$V2[idx],info$V4[idx],info$V5[idx]))
}


enrichout <- as.data.frame(enrichout)
for (idx in 2:9){
  enrichout[,idx] <- as.numeric(enrichout[,idx])
}
colnames(enrichout) <- c("classID","Motif_num_in_class","fg_target","fg_total","bg_target","bg_total","fg_frac","bg_frac","pval")
enrichout$pval[enrichout$fg_frac < enrichout$bg_frac] <- 1
write.table(enrichout,"./hg19_LoVo.class.txt",
            sep = "\t",row.names = F, col.names = T, quote = F)

enrichout2 <- as.data.frame(enrichout2)
for (idx in 2:8){
  enrichout2[,idx] <- as.numeric(enrichout2[,idx])
}
colnames(enrichout2) <- c("motifID","fg_target","fg_total","bg_target","bg_total","fg_frac","bg_frac","pval",
                          "motif_name","motif_class","motif_family")
enrichout2$pval[enrichout2$fg_frac < enrichout2$bg_frac] <- 1
write.table(enrichout2,"./hg19_LoVo.eachmotif.txt",
            sep = "\t",row.names = F, col.names = T, quote = F)