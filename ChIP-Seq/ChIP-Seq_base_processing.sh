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

## Trimming chip/input: Trimmomatic
java -jar trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_R1.fastq ./trim/${sample_prefix}_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15
java -jar trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_input_R1.fastq ./trim/${sample_prefix}_input_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15

## Mapping chip/input: STAR
STAR --genomeDir /mnt/flatfiles/organisms/new_organism/mus_musculus/101/index_star --runThreadN 16 --readFilesIn ./trim/${sample_prefix}_R1_trim.fastq --outReadsUnmapped Fastx --outFileNamePrefix ./star/${sample_prefix} --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20 --alignEndsProtrude 10 ConcordantPair --alignMatesGapMax 2000 --limitOutSAMoneReadBytes 10000000 --outMultimapperOrder Random --sjdbOverhang 100 --alignEndsType EndToEnd --alignIntronMax 1 --alignSJDBoverhangMin 999 --outFilterMultimapNmax 1
STAR --genomeDir /mnt/flatfiles/organisms/new_organism/mus_musculus/101/index_star --runThreadN 16 --readFilesIn ./trim/${sample_prefix}_input_R1_trim.fastq --outReadsUnmapped Fastx --outFileNamePrefix ./star/${sample_prefix}_input --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20 --alignEndsProtrude 10 ConcordantPair --alignMatesGapMax 2000 --limitOutSAMoneReadBytes 10000000 --outMultimapperOrder Random --sjdbOverhang 100 --alignEndsType EndToEnd --alignIntronMax 1 --alignSJDBoverhangMin 999 --outFilterMultimapNmax 1

## Generate BigWig normalized to mapped reads per sample chip/input: python package deeptools 3.5.1 / bamCoverage
bamCoverage -b ./star/${sample_prefix}.bam -o ./star/${sample_prefix}.bw -p 16 --binSize 25 --smoothLength 75 --normalizeUsing RPKM --outFileFormat bigwig
bamCoverage -b ./star/${sample_prefix}_input.bam -o ./star/${sample_prefix}_input.bw -p 16 --binSize 25 --smoothLength 75 --normalizeUsing RPKM --outFileFormat bigwig



TODO:
all steps: with input
peak calling
recount
...



