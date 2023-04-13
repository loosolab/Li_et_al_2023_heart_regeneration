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

## Trimming: Trimmomatic
java -jar /mnt/software/x86_64/packages/trimmomatic/0.39/trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_R1.fastq ./trim/${sample_prefix}_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15
java -jar /mnt/software/x86_64/packages/trimmomatic/0.39/trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_input_R1.fastq ./trim/${sample_prefix}_input_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15






TODO:
all steps: with input
peak calling
recount
...











## Mapping: STAR
STAR --genomeDir /mnt/flatfiles/organisms/new_organism/mus_musculus/101/index_star --runThreadN 16 --readFilesIn ./trim/${sample_prefix}_R1_trim.fastq --outReadsUnmapped Fastx --outFileNamePrefix ./star/${sample_prefix} --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20 --alignEndsProtrude 10 ConcordantPair --alignMatesGapMax 2000 --limitOutSAMoneReadBytes 10000000 --outMultimapperOrder Random --sjdbOverhang 100 --alignEndsType EndToEnd --alignIntronMax 1 --alignSJDBoverhangMin 999 --outFilterMultimapNmax 1

## Deduplication: PICARD
java -Xmx4g -jar picard.jar MarkDuplicates I=./star/${sample_prefix}.bam O=./star/${sample_prefix}_nodup.bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

## Identify ribosomal reads: samtools + bbtools + silva db
samtools bam2fq ./star/${sample_prefix}_nodup.bam >./star/${sample_prefix}_nodup.bam.fastq
bbduk.sh in=./star/${sample_prefix}_nodup.bam.fastq out=./star/${sample_prefix}_nodup_norrna.bam.fastq outm=./star/${sample_prefix}_nodup_rrna.bam.fastq ref=/mnt/flatfiles/rrna_silva/ribokmers.fa k=31 hdist=0 stats=./star/${sample_prefix}_rrnacheck_stats_detail.txt -eoom threads=16
cat ./star/${sample_prefix}_nodup_norrna.bam.fastq | sed -n '1~4s/^@/>/p;2~4p' | grep ">" | sed 's/>//g' | sed 's/\/1$//g' | sed 's/\/2$//g' | awk '!a[$0]++' >./star/${sample_prefix}_nodup.bam.rrna.txt

## Identify chrM ids: samtools
echo -e "chrM" | xargs --max-procs=8 samtools view ./star/${sample_prefix}_nodup.bam | cut -f1 | awk '!a[$0]++' >./star/${sample_prefix}_nodup.bam.chrm.txt

## Remove rRNA and chrM reads: PICARD
cat ./star/${sample_prefix}_nodup.bam.rrna.txt >./star/${sample_prefix}_nodup.bam.remove.txt
cat ./star/${sample_prefix}_nodup.bam.chrm.txt >>./star/${sample_prefix}_nodup.bam.remove.txt
java -Xmx24g -jar picard.jar FilterSamReads I=./star/${sample_prefix}_nodup.bam O=./star/${sample_prefix}_nodup_filter.bam READ_LIST_FILE=./star/${sample_prefix}_nodup.bam.remove.txt FILTER=excludeReadList

## Generate BigWig normalized to mapped reads per sample: python package deeptools 3.5.1 / bamCoverage
bamCoverage -b ./star/${sample_prefix}_nodup_filter.bam -o ./star/${sample_prefix}_nodup_filter.bw -p 16 --binSize 25 --smoothLength 75 --normalizeUsing RPKM --outFileFormat bigwig
