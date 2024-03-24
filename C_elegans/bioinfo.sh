######################################
####   XR-seq
######################################

# Cutting the adapter with Cutadapt
for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-cutadapt.out" --wrap="cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -m 1 -o ${SAMPLE}_trimmed.fastq ${SAMPLE_DIR}/${SAMPLE}.fastq.gz"
done

# Merging identical reads
for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --mem=8g -t 120 --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-fastx.out" --wrap="fastx_collapser -v -i ${SAMPLE}_trimmed.fastq -o ${SAMPLE}_trimmed.fasta -Q33"
done


# Genome alignment
for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bowtie2.out" --wrap="bowtie2 -x ${BOWTIE2_IND} -f ${SAMPLE}_trimmed.fasta -S ${SAMPLE}_trimmed.sam --very-sensitive"
done

# Convert to sorted BAM
for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --mem=16g --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-samtools_sort.out" --wrap="samtools sort -o ${SAMPLE}_trimmed_sorted.bam ${SAMPLE}_trimmed.sam"
done

# Convert to BED
for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-bamtobed.out" --wrap="bedtools bamtobed -i ${SAMPLE}_trimmed_sorted.bam > ${SAMPLE}_trimmed_sorted.bed"
done

# Filter by length
if [[ -n "$MIN" ]] && [[ -n "$MAX" ]]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-length_filter.out" --wrap="awk '{if(\$3-\$2>=$MIN && \$3-\$2<=$MAX){print}}' ${SAMPLE}_trimmed_sorted.bed > ${SAMPLE}_filtered.bed"
  done

######################################
####   RNA-seq
######################################

# STAR alignment (paired-end RNA sequencing), sam to bam, filter, and sort
STAR --genomeDir elegans_STAR --readFilesIn wtL1R1_1.fq wtL1R1_2.fq --outFilterIntronMotifs RemoveNoncanonicalUnannotated –runThreadN 2 --outFileNamePrefix wt_
samtools view -bS wt_Aligned.out.sam > wt_Aligned.out.bam
samtools sort -O bam -o wt_Aligned.out.sorted.bam wt_Aligned.out.bam
samtools index wt_Aligned.out.sorted.bam

STAR --genomeDir elegans_STAR --readFilesIn xpcL1R1_1.fq xpcL1R1_2.fq --outFilterIntronMotifs RemoveNoncanonicalUnannotated –runThreadN 2 --outFileNamePrefix xpc_
samtools view -bS xpc_Aligned.out.sam > xpc_Aligned.out.bam
samtools sort -O bam -o xpc_Aligned.out.sorted.bam xpc_Aligned.out.bam
samtools index xpc_Aligned.out.sorted.bam
