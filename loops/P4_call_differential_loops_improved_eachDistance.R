#---------------------------------------------------------------------------
#predefine parameters
res <- 10000
#---------------------------------------------------------------------------
#load raw counts
samples <- c("/siCtrl.1_10000.loop.counts.txt",
             "siCtrl.2_10000.loop.counts.txt",
             "siNIPBL.1_10000.loop.counts.txt",
             "siNIPBL.2_10000.loop.counts.txt")
counts1 <- NULL
for (i in 1:length(samples)){
  count.tmp <- read.table(samples[i],sep = "\t",header = F)
  counts1 <- cbind(counts1,count.tmp[,8])
}
raw.counts <- cbind(count.tmp[,1:6],counts1)
colnames(raw.counts)[7:10] <- c("siCtrl1","siCtrl2","siNIPBL1","siNIPBL2")


#---------------------------------------------------------------------------
#load the total counts at each certain distance to correct the bias of power-law difference across samples
pw <- c("siCtrl.1_10000.totalcount.at.a.distance",
        "siCtrl.2_10000.totalcount.at.a.distance",
        "siNIPBL.1_10000.totalcount.at.a.distance",
        "siNIPBL.2_10000.totalcount.at.a.distance")
pwlaw <- NULL
for (i in 1:length(pw)){
  pw.file <- read.table(paste0("./",pw[i]),sep = "\t",header = F)
  pwlaw <- cbind(pwlaw,pw.file[,2])
}


#---------------------------------------------------------------------------
#calculate the pvals with binomial distribution
pval.binomial <- NULL
for (i in 1:nrow(raw.counts)){
  bin_pos <- (raw.counts[i,5]-raw.counts[i,2])%/%res + 1
  p <- sum(pwlaw[bin_pos,1:2])/sum(pwlaw[bin_pos,1:4])
  total.count <- sum(raw.counts[i,7:10])
  query.count <- sum(raw.counts[i,7:8])
  scaled.val <- pwlaw[bin_pos,1:4]/(prod(pwlaw[bin_pos,1:4])^(1/4))
  norm.count <- unname(raw.counts[i,7:10]/scaled.val)
  if (total.count == 0){
    tmp.p <- 1
  } else {
    tmp.p <- min(pbinom(query.count,total.count,p),pbinom(query.count-1,total.count,p,lower.tail = F))
  }
  #-------------------------------------------------------------------------
  #background pvals
  p.bg1 <- sum(pwlaw[bin_pos,c(1,3)])/sum(pwlaw[bin_pos,1:4])
  query.count.bg1 <- sum(raw.counts[i,c(7,9)])
  if (total.count == 0){
    tmp.p.bg1 <- 1
  } else {
    tmp.p.bg1 <- min(pbinom(query.count.bg1,total.count,p.bg1),pbinom(query.count.bg1-1,total.count,p.bg1,lower.tail = F))
  }  
  p.bg2 <- sum(pwlaw[bin_pos,c(1,4)])/sum(pwlaw[bin_pos,1:4])
  query.count.bg2 <- sum(raw.counts[i,c(7,10)])
  if (total.count == 0){
    tmp.p.bg2 <- 1
  } else {
    tmp.p.bg2 <- min(pbinom(query.count.bg2,total.count,p.bg2),pbinom(query.count.bg2-1,total.count,p.bg2,lower.tail = F))
  } 
  pval.binomial <- rbind(pval.binomial, c(norm.count,tmp.p,tmp.p.bg1,tmp.p.bg2))
}

pvals <- as.data.frame(pval.binomial)
colnames(pvals) <- c("norm.siCtrl1","norm.siCtrl2","norm.siNIPBL1","norm.siNIPBL2","true.pval","bg.pval1","bg.pval2")

for (i in 1:ncol(pvals)) {
  pvals[,i] <- as.numeric(pvals[,i])
}

integrated.result <- cbind(raw.counts,pvals)
integrated.result$total.count <- apply(integrated.result[,7:10],1,sum)
integrated.result$logfc <- log2(apply(integrated.result[,13:14],1,sum)/apply(integrated.result[,11:12],1,sum))
hist(integrated.result$total.count,breaks = 50)


library(ggplot2)
ggplot(integrated.result,aes(y=-log10(true.pval),x = total.count)) +
  geom_point() + geom_rug() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

ggplot(integrated.result,aes(y=-log10(bg.pval2),x = total.count)) +
  geom_point() + geom_rug() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

ggplot(integrated.result,aes(y=-log10(bg.pval1),x = total.count)) +
  geom_point() + geom_rug() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

#investigate length of loops
ggplot(integrated.result,aes(y=-log10(true.pval),x = logfc,color = V5-V2)) +
  geom_point() + geom_rug(col = "#FF8000",alpha=0.1) + scale_colour_gradient(limits = c(0,750000)) + xlim(-2,2)
#investigate differences
integrated.result$label <- "Unchanged"

FDR <- 0.05
##determine the p value threshold according to FDR
p.thresh.vector <- rev(seq(0.001,0.05,by = 0.001))
for (x in p.thresh.vector){
  ture.fdr1 <-  sum(x > integrated.result$bg.pval1)/(sum(x > integrated.result$bg.pval1) + sum(x > integrated.result$true.pval))
  ture.fdr2 <-  sum(x > integrated.result$bg.pval2)/(sum(x > integrated.result$bg.pval2) + sum(x > integrated.result$true.pval))
  if (ture.fdr1 <= FDR & ture.fdr2 <= FDR){
    pval_thres <- x
    break
  } else{
    pval_thres <- x
  }
}

integrated.result$label[integrated.result$logfc > log2(1.5) & integrated.result$true.pval < pval_thres ] <- "Increased"
integrated.result$label[integrated.result$logfc < -log2(1.5) & integrated.result$true.pval < pval_thres] <- "Decreased"
integrated.result$label <- factor(integrated.result$label,levels = c("Increased","Decreased","Unchanged"))

write.table(integrated.result, "./diff_loops_siCtrlVSsiNIPBL_improved.txt",sep = "\t",row.names = F,
          quote = F, col.names = T)

pdf("./siNIPBL_diffloop.pdf",width = 7,height = 5)
ggplot(integrated.result,aes(y=-log10(true.pval),x = logfc,colour = label)) +
  geom_point() +xlim(-3,3) + scale_color_manual(values=c("#FF8000","#6A5ACD","#C0C0C0")) +
  annotate("text",x = -1.5,y=12,label = paste0("Num:",as.character(sum(integrated.result$label == "Decreased"))),colour = "#6A5ACD")+
  annotate("text",x = 1.5,y=12,label = paste0("Num:",as.character(sum(integrated.result$label == "Increased"))),colour = "#FF8000")+
  annotate("text",x = 0,y=2,label = paste0("Num:",as.character(sum(integrated.result$label == "Unchanged"))),colour = "black")
dev.off()
#geom_rug(col = "#A39480",alpha=0.1)+xlim(-3.5,3.5)+

#build a new dataframe
new.df <- data.frame(pvals = c(integrated.result$true.pval,integrated.result$bg.pval1,integrated.result$bg.pval2))
new.df$label <- "True.pval"
new.df$label[(nrow(integrated.result)+1):(nrow(integrated.result)*2)] <- "Bg1.pval"
new.df$label[(nrow(integrated.result)*2+1):(nrow(integrated.result)*3)] <- "Bg2.pval"
new.df$label <- factor(new.df$label,levels = c("True.pval","Bg1.pval","Bg2.pval"))
ggplot(new.df,aes(x=pvals,fill = label)) + geom_histogram(alpha = 0.5,position="identity",breaks = seq(0,0.5,0.01))


#
integrated.result$V1 <- factor(integrated.result$V1)
x1 <- table(integrated.result$V1[integrated.result$label == "Decreased"])
x2 <- table(integrated.result$V1[integrated.result$label == "Increased"])
x3 <- table(integrated.result$V1[integrated.result$label == "Unchanged"])
loopsinchrom <- cbind(x1,x2,x3)
loopsinchrom2 <- as.data.frame(loopsinchrom[1:22,])
loopsinchrom2$name <- rownames(loopsinchrom2)
loopsinchrom2$up.ratio <- loopsinchrom2$x1/loopsinchrom2$x3
loopsinchrom2$down.ratio <- loopsinchrom2$x2/loopsinchrom2$x3
ggplot(loopsinchrom2,aes(x=up.ratio,y=down.ratio,label = name)) + geom_point() + geom_text(nudge_x = 0.003,nudge_y = 0.001)+
  xlim(0,0.25) + ylim(0,0.035)
