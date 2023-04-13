ibrary(DESeq2)
options(bitmapType='cairo')

data = read.table("counts.matrix", header=T, row.names=1, com='', sep="\t", quote="", check.names=F)
col_ordering = c(1,2,3,4)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=5,]
cat("\nmatrix:\n")
colnames(rnaseqMatrix)
conditions = data.frame(conditions=factor(c("ctrl","ctrl","ko","ko")))
rownames(conditions) = colnames(rnaseqMatrix)
cat("\ndesign:\n")
conditions
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rnaseqMatrix, colData = conditions, design = ~ conditions)

dds = DESeq(ddsFullCountTable)
res = results(dds) 

baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "ctrl"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "ko"])

res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
		
##adds an id column headline for column 0
res = cbind(id=rownames(res), as.data.frame(res))

##replace "NA" values in padj column by 1
res$padj[is.na(res$padj)] <- 1
	
## output results; set row.names to false to accomodate change above
write.table(as.data.frame(res[order(res$pvalue),]), file='counts.matrix.ctrl_vs_ko-tet3.results', sep='	', quote=FALSE, row.names=F)
