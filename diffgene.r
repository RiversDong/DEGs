library(DESeq2)
exp_matrix_f <- "/home/chuand/new_gene/result/summary"
bampath <- "X.home.chuand.new_gene.data.mergebam."
args <- commandArgs(T)
metafile <- args[1]
differgene <- args[2]
allfoldchange <- args[3]

exp_matrix <- read.table(exp_matrix_f, header=TRUE, skip=1, row.names=1, sep="\t")
colnames(exp_matrix) <- gsub(bampath, "", colnames(exp_matrix), fixed=TRUE)
colnames(exp_matrix) <- gsub(".bam", "", colnames(exp_matrix), fixed=TRUE)
metadata <- read.table(metafile, header=TRUE, sep = "\t")
ids <- metadata$sampleId
exp_matrix <- exp_matrix[, ids]
countdata <- exp_matrix[,1:ncol(exp_matrix)]
rownames(metadata) <- metadata$sampleId

# construct the dds object
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design=~group)
dds <- DESeq(dds) # standardization
res <- results(dds, pAdjustMethod="fdr")
diff_gene <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(x=as.data.frame(diff_gene), file=differgene, sep="\t", quote=F, col.names=NA)
write.table(x=as.data.frame(res), file=allfoldchange, sep="\t", quote=F, col.names=NA)
