# STAR alignment

while read fastq
do
sbatch -t 0-23:00 --mem 60000 --wrap="/nas/longleaf/home/yuchaoj/yuchaojlab/program/STAR-2.5.3a/bin/Linux_x86_64/STAR --genomeDir /nas/longleaf/home/yuchaoj/yuchaojlab/lib --readFilesCommand zcat --readFilesIn "$fastq"_1.fq.gz "$fastq"_2.fq.gz --outFilterIntronMotifs RemoveNoncanonicalUnannotated –runThreadN 2 --outFileNamePrefix ~/scratch/RNAseq/"$fastq"_" --job-name $fastq
done < fastq.list

# 2) sam to bam, filter, sort and build index
while read fastq
do
sbatch -t 0-10:00 --mem 30000 --wrap="~/yuchaojlab/bin/samtools view -bS "$fastq"_Aligned.out.sam > "$fastq"_Aligned.out.bam; perl ~/yuchaojlab/bin/filter_sam_v2.pl "$fastq"_Aligned.out.bam "$fastq"_Aligned.out.filtered.sam; ~/yuchaojlab/bin/samtools view -bS "$fastq"_Aligned.out.filtered.sam > "$fastq"_Aligned.out.filtered.bam; java -Xmx30G -jar ~/yuchaojlab/bin/picard.jar SortSam INPUT="$fastq"_Aligned.out.filtered.bam OUTPUT="$fastq"_Aligned.out.filtered.sorted.bam SORT_ORDER=coordinate; java -jar ~/yuchaojlab/bin/picard.jar BuildBamIndex I="$fastq"_Aligned.out.filtered.sorted.bam"  --job-name $fastq
done < fastq.list



# below is using featurecounts to get read counts

while read fastq; do
sbatch -t 0-15:00 --mem 40000 --wrap="/nas/longleaf/home/yuchaoj/yuchaojlab/program/subread-1.6.0-source/bin/featureCounts -t exon -g gene_id -a /nas/longleaf/home/yuchaoj/yuchaojlab/lib/Homo_sapiens.GRCh37.75.gtf -o "$fastq".featurecounts.txt "$fastq"_Aligned.out.filtered.sorted.bam"  --job-name $fastq
done < fastq.list

