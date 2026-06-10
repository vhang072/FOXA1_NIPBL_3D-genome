#load output from python script S2_call_Diff
data <- read.table("../dTAG_0h_vs_24h_20kb.dynamics",sep = "\t",
                   header = T)


new.df <- data.frame(pvals = c(data$TrueP,data$BgP1,data$BgP2))
new.df$label <- "True.pval"
new.df$label[(nrow(data)+1):(nrow(data)*2)] <- "Bg1.pval"
new.df$label[(nrow(data)*2+1):(nrow(data)*3)] <- "Bg2.pval"
new.df$label <- factor(new.df$label,levels = c("True.pval","Bg1.pval","Bg2.pval"))
new.df$pvals <- 10^(-1*abs(new.df$pvals))
library(ggplot2)
ggplot(new.df,aes(x=pvals,fill = label)) + geom_histogram(alpha = 0.5,position="identity",breaks = seq(0,0.5,0.01))

data$totalC <- apply(data[,4:7],1,sum)
data$log2fold <- log2((data$scaleS3 + data$scaleS4+1)/(data$scaleS1 + data$scaleS2+1))


FDR <- 0.05
##determine the p value threshold according to FDR
logp.thresh.vector <- seq(1.3, 12, by = 0.1)
for (x in logp.thresh.vector) {
    n_true <- sum(abs(data$TrueP) > x)
    n_bg1  <- sum(abs(data$BgP1) > x)
    n_bg2  <- sum(abs(data$BgP2) > x)

    fdr1 <- n_bg1 / max(n_bg1 + n_true, 1)
    fdr2 <- n_bg2 / max(n_bg2 + n_true, 1)

    if (fdr1 <= FDR & fdr2 <= FDR) {
        logp.thresh <- x
        break
    }
}

fc.thresh <- log2(1.5)
data$label <- "Unchanged"
data$label[(abs(data$TrueP)>logp.thresh) & (data$log2fold>fc.thresh) & (abs(data$BgP1)< -log10(0.05)) &
              (abs(data$BgP2)< -log10(0.05))] <- "Increased"
data$label[(abs(data$TrueP)>logp.thresh) & (data$log2fold< -1*fc.thresh) & (abs(data$BgP1)< -log10(0.05)) &
             (abs(data$BgP2)< -log10(0.05))] <- "Decreased"
data$label <- factor(data$label)
table(data$label)
diff.data <- data[data$label != "Unchanged",]
write.table(diff.data ,file="../dTAG_0h_vs_24h_20kb.diff_dynamics",
            sep = "\t",col.names = T, row.names = F, quote = F)

#---------------------------------------------------------------------------------------------------------
#QC
pvaldis <- NULL
cared_quantile <- c(0.5,0.7,0.9,0.93,0.96,0.99)
for (i in seq(200,2000,200)){
  f <- quantile(abs(data$BgP2)[(data$totalC >= i-200) & (data$totalC < i)],cared_quantile)
  pvaldis <- rbind(pvaldis,c(unname(f),sum((data$totalC >= i-200) & (data$totalC < i))))
}
#---------------------------------------------------------------------------------------------------------
ggplot(data,aes(x=totalC,y=abs(TrueP))) + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  xlim(0,1000) 

