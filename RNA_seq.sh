# Generate genome (read length 50bp)
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir Mus_STAR --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm38.91.gtf --sjdbOverhang 49

# STAR alignment (paired-end RNA sequencing)
STAR --genomeDir Mus_STAR --readFilesIn sample_R1.fastq sample_R2.fastq --outFilterIntronMotifs RemoveNoncanonicalUnannotated –runThreadN 2 --outFileNamePrefix sample_


# Sam to bam, filter, and sort
# filter_sam_v2.pl available in this GitHub page
samtools view -bS sample_Aligned.out.sam > sample_Aligned.out.bam
perl filter_sam_v2.pl sample_Aligned.out.bam sample_Aligned.out.filtered.sam
samtools view -bS sample_Aligned.out.filtered.sam > sample_Aligned.out.filtered.bam
java -Xmx30G -jar picard.jar SortSam INPUT=sample_Aligned.out.filtered.bam OUTPUT=sample_Aligned.out.filtered.sorted.bam SORT_ORDER=coordinate

# Build index
java -jar picard.jar BuildBamIndex I=sample_Aligned.out.filtered.sorted.bam

# Use featurecounts to get read counts across all genes
featureCounts -t exon -g gene_id -a Mus_musculus.GRCm38.91.gtf -o sample.featurecounts.txt sample_Aligned.out.filtered.sorted.bam" 

# Use SALMON to get transcript per million matrix (assembly free)
salmon quant -i transcripts_index -l A -1 sample_R1.fastq -2 sample_R2.fastq -o sample
