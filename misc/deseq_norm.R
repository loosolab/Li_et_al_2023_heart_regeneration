#!/usr/bin/env Rscript

#normalize a raw count matrix using deseq:
#DESeq: This normalization method is based on the hypothesis that most genes are not DE. A DESeq scaling factor for a
#given lane is computed as the median of the ratio, for each gene, of its read count over its geometric mean
#across all lanes. The underlying idea is that non-DE genes should have similar read counts across samples,
#leading to a ratio of 1. Assuming most genes are not DE, the median of this ratio for the lane provides an
#estimate of the correction factor that should be applied to all read counts of this lane to fulfill the hypothesis.
#By calling the estimateSizeFactors() and sizeFactors() functions in the DESeq Bioconductor package, this factor is computed for
#each lane, and raw read counts are divided by the factor associated with their sequencing lane.

#ATTENTION:
#type=ratio:
#	only rows without a 0 in all columns are used; outlier sample with (nearly) all 0 -> size factor will be set to NA
#	this is the deseq default, but i use poscounts as default here to avoid outlier/extreme sample issues

#format matrix:
# 1. line = headline
# 1. column = gene id
# rest: raw counts

#optional:
#format reps: condition_id sample_id (sample_id must be equals to the column names in matrix)
#m1	m1_1
#m1	m1_2
#m2	m2_1
#m2	m2_2

###############################################################################################################
# PARSE PARAMETERS

require('getopt');

##0=no, 1=required, 2=optional
##logical,integer,double,complex,character

options = matrix(c(
  'matrix','m',1,'character','Tab-delimited input matrix bearing count columns (1. row = column ids, 1. column = row ids).',
  'reps','r',1,'character','Optional: tab-delimited input file bearing association of conditions and samples (no headline, condition[TAB]sample). This is only used for plotting (to assign colors based on condition).',
  'outfile','o',2,'character','Output matrix (default=<input>.norm).',
  'type','t',2,'character','Type of normalization: ratio/poscounts. ratio = median ratio excluding features with 0 in at least one sample; poscounts = median ratio including all features; only works with R>=4.0 + DEseq>=1.30.0. Default: ratio',
  'dropouts','d',2,'double','Maximum ratio of dropouts when using normalization type ratio (-t ratio). Default: 0.6.',
  'help','h',0,'logical','Provides command line help.'),ncol=5,byrow=T)

opt = getopt(options);

# help was asked for.
if ( !is.null(opt$help) ) {
  cat(getopt(options, usage=TRUE),file = stderr());
  q(status=1);
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
if ( is.null(opt$matrix) ) {
  cat(getopt(options, usage=TRUE),file = stderr());
  q(status=2);
}

if ( is.null(opt$outfile) ) { 
	#fileprefix = sub("^(.*)[.].*", "\\1", opt$matrix)				#remove suffix from filename (everything after last ".") = output file prefix
	opt$outfile = paste(opt$matrix, ".norm", sep="")
}

if ( is.null(opt$type) ) {
  opt$type = "ratio"
}
if ( is.null(opt$dropouts) ) {								#maximum ratio of rows with >= 1 count of 0 in a column (=dropouts that will be ignored when using type=ratio)
  opt$dropouts = 0.6
}

opt$outfile_plot_raw=paste(opt$matrix, "_boxplot.pdf", sep="")
opt$outfile_plot_norm=paste(opt$outfile, "_boxplot.pdf", sep="")

##############################################################################################################

#options(stringsAsFactors = FALSE)	#to deal with potential issues in R 4 which sets this to true by default

library(DESeq2)
setwd(".")

##read in your data to make counts table
countsTable <- read.delim(opt$matrix, header=TRUE, stringsAsFactors=TRUE, row.names=1, check.names=F)
countsTable=round(countsTable)
countsTable[is.na(countsTable)] <- 0

rows_total=nrow(countsTable)
rows_dropout=sum(rowSums(countsTable == 0, na.rm=TRUE) > 0)		#nr. of rows with >= 1 count of 0 in a column
dropout_ratio=rows_dropout/rows_total
message("ratio of rows with a 0 in a column (=dropouts when using -t ratio): ",dropout_ratio)
if ( opt$type=="ratio" && dropout_ratio>opt$dropouts ) {
	message("nr. of rows: ", rows_total)
	message("nr. of rows with a 0 in a column: ",rows_dropout)
	message("ratio of rows with a 0 in a column (=dropouts): ",dropout_ratio)
	message("Error. ratio of rows with a 0 in a column > ",opt$dropouts,". please switch to type poscounts (-t poscounts). exit.")
	message("if this error happened inside a pipeline, please set deseqnorm_type=\"poscounts\" to use all features to compute factors OR disable the threshold with max_dropouts=\"1\" to only use features without a 0")
	stop("")
}

cat("\ncolumns:\n")
(condition <- factor(colnames(countsTable)))

#cat("\ncoldata\n")
coldata <- data.frame(row.names=colnames(countsTable), condition)

#cat("\ndds\n")
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countsTable), colData = coldata, design=~condition)

if (utils::packageVersion("DESeq2") < "1.30.0") {
	if ( opt$type=="poscounts" ) {
		stop("-t poscounts not available for DESeq2 < 1.30.0. please update package or switch to -t ratio. exit.")
	}
	dds <- estimateSizeFactors(dds)				#R 3 = type: ratio
} else {
	dds <- estimateSizeFactors(dds, type=opt$type)		#R 4 = type: ratio or poscounts
}    
factors=sizeFactors(dds)
factors

if (any(is.na(factors))) {	
	stop("deseqnorm failed to identify all size factors! there were no genes which are not 0 in all samples (outlier samples?) -> crash!\n")
}

scaledcounts = as.data.frame(counts(dds, normalized=TRUE))
scaledcounts = round(scaledcounts)
#cat("\nscaledcounts\n")
#head(scaledcounts)

##output
id<-data.frame(rownames(scaledcounts))		# making a data frame of the rownames from the scaledcounts data frame
colnames(id)<-c("id")				# renaming the column header to "id"
scaledcountsID<-cbind(id,scaledcounts)		# combining the columns (column binding)
#cat("\nscaledcountsID\n")
#head(scaledcountsID)
write.table(scaledcountsID, file=opt$outfile, sep="\t", quote=FALSE,  row.names=F)
