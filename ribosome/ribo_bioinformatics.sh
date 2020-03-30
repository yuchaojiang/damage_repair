# 1) cut adaptor sequence, remove reads too short (<= 10), discard untrimmed reads

while read fastq
do
sbatch -t 1-23:00 --mem 80000 --wrap="module load cutadapt; cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -m 22 -M 30 -o "$fastq".cu.fastq "$fastq".fastq" --job-name $fastq
done < fastq.list

# 2) index reference by BWA
srun -t 2-23:00 --mem 50000 --pty /bin/bash
~/yuchaojlab/bin/bwa index ribo.fa

# 3) alignment using BWA
# BWA manual page:  http://bio-bwa.sourceforge.net/bwa.shtml
while read fastq
do
sbatch -t 1-23:00 --mem 80000 --wrap="~/yuchaojlab/bin/bwa aln -t 16 ribo.fa "$fastq".cu.fastq > "$fastq".ribo.cu.sai; ~/yuchaojlab/bin/bwa samse ribo.fa "$fastq".ribo.cu.sai "$fastq".cu.fastq > "$fastq".ribo.cu.sam; ~/yuchaojlab/bin/samtools view -bS "$fastq".ribo.cu.sam > "$fastq".ribo.cu.bam" --job-name $fastq
done < fastq.list


# 4) filter, sort, dedup, index
# samtools filtering: -b output bam format
# -F 4: only output mapped reads
# -q 1: skip reads with mapq smaller than 1
while read fastq
do
sbatch -t 1-23:00 --mem 80000 --wrap="~/yuchaojlab/bin/samtools view -bF 4 -q 1 "$fastq".ribo.cu.bam > "$fastq".ribo.cu.filtered.bam; java -jar ~/yuchaojlab/bin/picard.jar SortSam I="$fastq".ribo.cu.filtered.bam O="$fastq".ribo.cu.filtered.sorted.bam SORT_ORDER=coordinate; java -jar ~/yuchaojlab/bin/picard.jar MarkDuplicates I="$fastq".ribo.cu.filtered.sorted.bam O="$fastq".ribo.cu.filtered.sorted.dedup.bam  M="$fastq".ribo.cu.filtered.sorted.dedup.metrics.txt; java -jar ~/yuchaojlab/bin/picard.jar BuildBamIndex I="$fastq".ribo.cu.filtered.sorted.dedup.bam" --job-name $fastq
done < fastq.list


# 5) no mismatch allowed, no gaps allowed
while read fastq
do
sbatch -t 1-23:00 --mem 80000 --wrap="module load bamtools; bamtools filter -tag XM:0 -in "$fastq".ribo.cu.filtered.sorted.dedup.bam -out "$fastq".ribo.cu.filtered.sorted.dedup.noMismatch.bam; bamtools filter -tag XO:0 -in "$fastq".ribo.cu.filtered.sorted.dedup.noMismatch.bam -out "$fastq".ribo.cu.filtered.sorted.dedup.noMismatch.noGap.bam; bamtools filter -tag X1:0 -in "$fastq".ribo.cu.filtered.sorted.dedup.noMismatch.noGap.bam -out "$fastq".ribo.cu.filtered.sorted.dedup.noMismatch.noGap.noSub.bam; java -jar ~/yuchaojlab/bin/picard.jar BuildBamIndex I="$fastq".ribo.cu.filtered.sorted.dedup.noMismatch.noGap.noSub.bam" --job-name $fastq
done < fastq.list



# 6) get total number of reads from fastq
srun -t 2-23:00 --mem 50000 --pty /bin/bash
while read fastq
do
echo $(cat $fastq.fastq | wc -l)/4|bc
done < fastq.list

while read fastq
do
echo $(cat $fastq.cu.fastq | wc -l)/4|bc
done < fastq.list

# get total number of mapped reads
while read fastq
do
~/yuchaojlab/bin/samtools view -F 0x904 -c $fastq.ribo.cu.filtered.bam
done < fastq.list

# get total number of mapped reads with no mismatches and no gaps
while read fastq
do
~/yuchaojlab/bin/samtools view -F 0x904 -c $fastq.ribo.cu.filtered.sorted.dedup.noMismatch.noGap.noSub.bam
done < fastq.list

