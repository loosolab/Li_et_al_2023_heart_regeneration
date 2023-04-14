## Reference GTF/FASTA: 
Ensembl mus musculus release 101 (mm10)

## Samples:
prefix	fastq_chip  fastq_input
ChIP-Cre_1	ChIP-Cre_1_R1.fastq	ChIP-Cre_1_input_R1.fastq
ChIP-Cre_2	ChIP-Cre_2_R1.fastq	ChIP-Cre_2_input_R1.fastq
ChIP-Cre_3	ChIP-Cre_3_R1.fastq	ChIP-Cre_3_input_R1.fastq
ChIP-KO_1	ChIP-KO_1_R1.fastq	ChIP-KO_1_input_R1.fastq
ChIP-KO_2	ChIP-KO_2_R1.fastq	ChIP-KO_2_input_R1.fastq
ChIP-KO_3	ChIP-KO_3_R1.fastq	ChIP-KO_3_input_R1.fastq

##
## Repeat the following steps per sample
##

## Trimming chip/input: Trimmomatic
java -jar trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_R1.fastq ./trim/${sample_prefix}_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15
java -jar trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_input_R1.fastq ./trim/${sample_prefix}_input_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15

## Mapping chip/input: STAR
STAR --genomeDir /mnt/flatfiles/organisms/new_organism/mus_musculus/101/index_star --runThreadN 16 --readFilesIn ./trim/${sample_prefix}_R1_trim.fastq --outReadsUnmapped Fastx --outFileNamePrefix ./star/${sample_prefix} --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20 --alignEndsProtrude 10 ConcordantPair --alignMatesGapMax 2000 --limitOutSAMoneReadBytes 10000000 --outMultimapperOrder Random --sjdbOverhang 100 --alignEndsType EndToEnd --alignIntronMax 1 --alignSJDBoverhangMin 999 --outFilterMultimapNmax 1
STAR --genomeDir /mnt/flatfiles/organisms/new_organism/mus_musculus/101/index_star --runThreadN 16 --readFilesIn ./trim/${sample_prefix}_input_R1_trim.fastq --outReadsUnmapped Fastx --outFileNamePrefix ./star/${sample_prefix}_input --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20 --alignEndsProtrude 10 ConcordantPair --alignMatesGapMax 2000 --limitOutSAMoneReadBytes 10000000 --outMultimapperOrder Random --sjdbOverhang 100 --alignEndsType EndToEnd --alignIntronMax 1 --alignSJDBoverhangMin 999 --outFilterMultimapNmax 1

## Deduplication chip/input: PICARD
java -Xmx4g -jar picard.jar MarkDuplicates I=./star/${sample_prefix}.bam O=./star/${sample_prefix}_nodup.bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
java -Xmx4g -jar picard.jar MarkDuplicates I=./star/${sample_prefix}_input.bam O=./star/${sample_prefix}_input_nodup.bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

## Generate BigWig normalized to mapped reads per sample chip/input: python package deeptools 3.5.1 / bamCoverage
bamCoverage -b ./star/${sample_prefix}_nodup.bam -o ./star/${sample_prefix}.bw -p 16 --binSize 25 --smoothLength 75 --normalizeUsing RPKM --outFileFormat bigwig
bamCoverage -b ./star/${sample_prefix}_input_nodup.bam -o ./star/${sample_prefix}_input.bw -p 16 --binSize 25 --smoothLength 75 --normalizeUsing RPKM --outFileFormat bigwig

## Call peaks: Music
samtools view ./star/${sample_prefix}_nodup.bam | MUSIC -preprocess SAM stdin ./music/${sample_prefix}/chip/
samtools view ./star/${sample_prefix}_input_nodup.bam | MUSIC -preprocess SAM stdin ./music/${sample_prefix}/input/
MUSIC -sort_reads ./music/${sample_prefix}/chip/ ./music/${sample_prefix}/chip/sorted
MUSIC -sort_reads ./music/${sample_prefix}/input/ ./music/${sample_prefix}/input/sorted
MUSIC -remove_duplicates ./music/${sample_prefix}/chip/sorted 2 ./music/${sample_prefix}/chip/dedup
MUSIC -remove_duplicates ./music/${sample_prefix}/input/sorted 2 ./music/${sample_prefix}/input/dedup
MUSIC -get_multiscale_punctate_ERs -chip ./music/${sample_prefix}/chip/dedup -control ./music/${sample_prefix}/input/dedup -l_p 1500 -mapp /mnt/flatfiles/organisms/mouse/mm10_GRCm38/music_mappability/100bp -l_mapp 100 -q_val 0.2 -sigma 0

##
## Combine per-sample results
##

## Create union peaks
cat ./music/*/*.bed | bedtools merge -header -d 50 >matrix/union.bed

## Recount union peaks for all samples: Run per sample
bigWigAverageOverBed ./star/${sample_prefix}.bw ./matrix/union.bed ./matrix/bed_recount/${sample_prefix}.union.counts.txt
cut -f7 ./matrix/bed_recount/${sample_prefix}.union.counts.txt >./matrix/bed_recount/merge/${sample_prefix}

## Merge union peak counts per sample to count matrix
paste ./matrix/bed_recount/merge/* >counts.matrix

## RNA differential analysis: R with DESeq2 1.26.0 package
## de_analysis.R script has to be adapted to include samples/conditions of the respective dataset
R --vanilla -q < de_analysis.R

## Normalize matrix for downstream plots: R with DESeq2 1.26.0 package, =counts.matrix.norm
deseq_norm.R -m counts.matrix


