# Trim 3' end adaptor sequence and discard short (<= 10bp and untrimmed reads
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -m 10 -o sample.cu.fastq sample.fastq

# Align using BWA (single-end sequencing, read length 50bp)
bwa aln -t 16 Mus_musculus.GRCm38.dna.primary_assembly.fa sample.cu.fastq > sample.cu.sai
bwa samse Mus_musculus.GRCm38.dna.primary_assembly.fa sample.cu.sai sample.cu.fastq > sample.cu.sam

# Sam to bam, then sort, dedup, index
samtools view -bS sample.cu.sam > sample.cu.bam
java -jar picard.jar SortSam I=sample.cu.bam O=sample.cu.sorted.bam SORT_ORDER=coordinate
java -jar picard.jar MarkDuplicates I=sample.cu.sorted.bam O=sample.cu.sorted.dedup.bam  M=sample.cu.sorted.dedup.metrics.txt
java -jar picard.jar BuildBamIndex I=sample.cu.sorted.dedup.bam
