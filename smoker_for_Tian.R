library(limma)
library(DESeq2)
#######


#set working directory
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/R/2022_for_Dominique")

#read in files 
raw.genecounts = read.table (file = "bulk_RNAseq_raw_counts.txt.minRib.txt.PC.txt",header = T,sep = "\t",row.names = 1)

#correct for UMI saturation the log = ln!!!! no log10 no log2 its natural log
raw.genecounts=round(-4096*(log(1-(raw.genecounts/4096))))

#check header
raw.genecounts[1:10,1:10]

#Scale before quantile normalization
raw.genecounts=t(t(raw.genecounts)/colSums(raw.genecounts))*100000

#quantile normalize
raw.genecounts=round(limma::normalizeQuantiles(raw.genecounts))


metadata <- read.table(file ="bulk_RNAseq_clinical_data_MM_selected.txt" ,header = T ,sep = "\t",row.names = 1) ## read metadata NOTE sep="\t"
head(metadata)
metadata = na.omit(metadata)   ## delete the na values 

raw.genecounts_selcted <- as.data.frame(raw.genecounts[ ,rownames(metadata)] )  ## select the gene from rawgenecounts

raw.genecounts_selcted=na.omit(raw.genecounts_selcted)   ## delete the na values 

dds <- DESeqDataSetFromMatrix(raw.genecounts_selcted ,metadata , design = ~ smokercurrent) 
dds <- DESeq(dds)

res <- results(dds)
head(res)
plotMA (dds, ylim=c(-1,1))   ## note : ylim

##  to see DEGs 

resOrdered <- res[order(res$padj),]  ## order the gene by padj
head(resOrdered)
resOrdered <- res[order(res$pvalue) ,] ## order the gene by pvalue
head(resOrdered)

sig <- res[!is.na(res$pvalue) & res$pvalue <0.05 &  abs(res$log2FoldChange > 0.4) ,]  ## select important different expression gene
head(sig)

##  write results to tab delimited  file
sig <- res[ !is.na(res$pvalue) & res$pvalue <0.05 &  abs(res$log2FoldChange ) > 1, ]

write.table(as.data.frame(sig), file= "current_smoke_results.xls" ,sep = "\t")    ## write results to excel file 



## load libraries for the heat map
library("gplots") 

selected <- rownames(sig)  ## selecte the gene 

# colors of the heat map
hmcol <- topo.colors(100)
hmcol <- heat.colors(100)

# heatmap
heatmap(    log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]+0.1)   , col = hmcol, scale="row" ) 


vignette("DESeq2")  ###  Analyzing RNA-seq data with DESeq2 from mike

plotDispEsts(dds)  ## Dispersion plot


## PCAplot
vsdata <- vst(dds,blind = F)
rld <- rlogTransformation(dds)
plotPCA(vsdata, intgroup = "smokercurrent") 

sum(res$pvalue < 0.1 ,na.rm = TRUE)  ## see how many gene pvalue<0.1 ,omit na

plotCounts(dds,gene = which.min(res$padj) ,intgroup = "smokercurrent")  ## plot counts of single gene
