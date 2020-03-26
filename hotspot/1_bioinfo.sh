# XR-seq pipeline
reference=~/yuchaojlab/lib/ucsc.hg19.fasta
while read fastq
do
sbatch -t 0-15:00 --mem 80000 --wrap="module load cutadapt; cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -m 22 -M 30 -o "$fastq".cu.fastq "$fastq".fastq; ~/yuchaojlab/bin/bwa aln -t 16 "$reference" "$fastq".cu.fastq > "$fastq".cu.sai; ~/yuchaojlab/bin/bwa samse "$reference" "$fastq".cu.sai "$fastq".cu.fastq > "$fastq".cu.sam; ~/yuchaojlab/bin/samtools view -bS "$fastq".cu.sam > "$fastq".cu.bam; ~/yuchaojlab/bin/samtools view -bF 4 -q 1 "$fastq".cu.bam > "$fastq".cu.filtered.bam; java -jar ~/yuchaojlab/bin/picard.jar SortSam I="$fastq".cu.filtered.bam O="$fastq".cu.filtered.sorted.bam SORT_ORDER=coordinate; java -jar ~/yuchaojlab/bin/picard.jar BuildBamIndex I="$fastq".cu.filtered.sorted.bam" --job-name $fastq
done < fastq.list
