## Reference GTF/FASTA: Ensembl mus musculus release 101 (mm10)
## place fastq files in subfolder ./raw/
## prepare STAR genome index as described in manual

## Example for samples (different per dataset):
prefix	fastq
ctrl_1	ctrl_1_R1.fastq
ctrl_2	ctrl_2_R1.fastq
ctrl_3	ctrl_3_R1.fastq
ko_1	ko_1_R1.fastq
ko_2	ko_2_R1.fastq
ko_3	ko_3_R1.fastq

##
## Repeat the following steps per sample
##

## Trimming: Trimmomatic
java -jar /mnt/software/x86_64/packages/trimmomatic/0.39/trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_R1.fastq ./trim/${sample_prefix}_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15

## Mapping: STAR
/mnt/software/x86_64/packages/star/2.7.7a/STAR --genomeDir /mnt/flatfiles/organisms/new_organism/mus_musculus/101/index_star --runThreadN 16 --readFilesIn ./trim/${sample_prefix}_R1_trim.fastq --outReadsUnmapped Fastx --outFileNamePrefix ./star/${sample_prefix} --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.9 --outFilterMatchNminOverLread 0.9 --outFilterMatchNmin 20 --alignIntronMax 200000 --alignMatesGapMax 2000 --alignEndsProtrude 10 ConcordantPair --outMultimapperOrder Random --limitOutSAMoneReadBytes 10000000 --sjdbOverhang 100 --outFilterMultimapNmax 1

## Deduplication: PICARD
java -Xmx4g -jar /mnt/software/x86_64/packages/picard/2.21.7/picard.jar MarkDuplicates I=./star/${sample_prefix}.bam O=./star/${sample_prefix}_nodup.bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

## Identify ribosomal reads: samtools + bbtools + silva db
/mnt/software/x86_64/packages/samtools/1.11/bin/samtools bam2fq ./star/${sample_prefix}_nodup.bam >./star/${sample_prefix}_nodup.bam.fastq
/mnt/software/x86_64/packages/bbtools/38.90/bbduk.sh in=./star/${sample_prefix}_nodup.bam.fastq out=./star/${sample_prefix}_nodup_norrna.bam.fastq outm=./star/${sample_prefix}_nodup_rrna.bam.fastq ref=/mnt/flatfiles/rrna_silva/ribokmers.fa k=31 hdist=0 stats=./star/${sample_prefix}_rrnacheck_stats_detail.txt -eoom threads=16
cat ./star/${sample_prefix}_nodup_norrna.bam.fastq | sed -n '1~4s/^@/>/p;2~4p' | grep ">" | sed 's/>//g' | sed 's/\/1$//g' | sed 's/\/2$//g' | awk '!a[$0]++' >./star/${sample_prefix}_nodup.bam.rrna.txt

## Identify chrM ids: samtools
echo -e "chrM" | xargs --max-procs=8 /mnt/software/x86_64/packages/samtools/1.11/bin/samtools view ./star/${sample_prefix}_nodup.bam | cut -f1 | awk '!a[$0]++' >./star/${sample_prefix}_nodup.bam.chrm.txt

## Remove rRNA and chrM reads: PICARD
cat ./star/${sample_prefix}_nodup.bam.rrna.txt >./star/${sample_prefix}_nodup.bam.remove.txt
cat ./star/${sample_prefix}_nodup.bam.chrm.txt >>./star/${sample_prefix}_nodup.bam.remove.txt
java -Xmx24g -jar /mnt/software/x86_64/packages/picard/2.21.7/picard.jar FilterSamReads I=./star/${sample_prefix}_nodup.bam O=./star/${sample_prefix}_nodup_filter.bam READ_LIST_FILE=./star/${sample_prefix}_nodup.bam.remove.txt FILTER=excludeReadList

## Generate BigWig normalized to mapped reads per sample: python package deeptools 3.5.1 / bamCoverage
/mnt/software/x86_64/packages/python/3.8.0-miniconda-4.9.2-buster/bin/bamCoverage -b ./star/${sample_prefix}_nodup_filter.bam -o ./star/${sample_prefix}_nodup_filter.bw -p 16 --binSize 25 --smoothLength 75 --normalizeUsing RPKM --outFileFormat bigwig

## Count reads per genes: subread featureCounts
/mnt/software/x86_64/packages/subread/1.6.5/bin/featureCounts -T 16 -o ./star/${sample_prefix}_fcounts.txt -t exon -g gene_id -s 1 -a mus_musculus.101.mainChr.gtf ./star/${sample_prefix}_nodup_filter.bam


##
## Combine per-sample results
##

## Join all sample counts to raw count matrix: =counts.matrix
merge_counts.pl ./star/*_fcounts.txt >counts.matrix

## RNA differential analysis: R with DESeq2 1.26.0 package
## de_analysis.R script has to be adapted to include samples/conditions of the respective dataset
R --vanilla -q < de_analysis.R

## Normalize matrix for downstream plots: R with DESeq2 1.26.0 package, =counts.matrix.norm
deseq_norm.R -m counts.matrix
