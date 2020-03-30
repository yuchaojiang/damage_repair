# Discard reads with adaptor sequence at 5' end, remove reads with length <= 10
cutadapt -g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT --discard-trimmed -m 10 -o sample.cu.fastq sample.fastq

# Index mouse genome by BWA
bwa index Mus_musculus.GRCm38.dna.primary_assembly.fa

# BWA alignment (single-end with 50bp read length)
bwa aln -t 16 Mus_musculus.GRCm38.dna.primary_assembly.fa sample.cu.fastq > sample.cu.sai 
bwa samse Mus_musculus.GRCm38.dna.primary_assembly.fa sample.cu.sai sample.cu.fastq > sample.cu.sam

# Sam to bam
samtools view -bS sample.cu.sam > sample.cu.bam

# Sort bam
java -jar picard.jar SortSam I=sample.cu.bam O=sample.cu.sorted.bam SORT_ORDER=coordinate

# Deduplication
java -jar picard.jar MarkDuplicates I=sample.cu.sorted.bam O=sample.cu.sorted.dedup.bam  M=sample.cu.sorted.dedup.metrics.txt

# Index
java -jar picard.jar BuildBamIndex I=sample.cu.sorted.dedup.bam
