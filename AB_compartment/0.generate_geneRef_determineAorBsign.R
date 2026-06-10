#choose only gene record lines from gene annotation gtf, for example gencode GTF
gene_anno <- read.table("./gencode.vxx.hg19.onlyGenes.txt",
                        sep = "\t",header = F)
library(stringr)
x <- matrix(0,25000,23)
res <- 10000

for (i in 1:nrow(gene_anno)){
  chrom = str_split(gene_anno[i,1],"r")[[1]][2]
  if (gene_anno[i,1] != "chrY" & gene_anno[i,1] != "chrM" & nchar(gene_anno[i,1])<=5 ){
    if (chrom == "X"){
      start = floor(gene_anno[i,4]/res)
      end = floor(gene_anno[i,5]/res)
      x[start:end,23] <- x[start:end,23] + 1/(end-start+1)
    } else {
      start = floor(gene_anno[i,4]/res)
      end = floor(gene_anno[i,5]/res)
      x[start:end,as.numeric(chrom)] <- x[start:end,as.numeric(chrom)] + 1/(end-start+1)
    }
  }
}
  
write.table(x,"./gencode.vxx.hg19_ABcompRef.txt",
            sep = "\t",row.names = F, col.names = F,quote = F)  