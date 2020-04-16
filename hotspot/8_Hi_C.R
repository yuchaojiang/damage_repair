#####################
### FIRE calling using FIREcaller ###
##############################
library("FIREcaller")

### Hi-C input file
# The Hi-C input file prefix.list is defined according to the naming convention of the NxN matrices, and the files for the contact matrices must have the naming convention "${prefix}_chr${number}.gz".

prefix.list <- c('IMR90_40Kb_raw_count_matrix')

### Define the genome build
gb<-'hg19'

### Define the name of the mappability file
map_file<-'F_GC_M_HIND3_40KB_HG19.txt.gz'

### Define whether to remove MHC region
# The MHC regions defined in FIREcaller are: hg19 chr6 28477797-33448354
rm_mhc <- TRUE

### Call FIREs and SuperFIREs
FIREcaller(prefix.list, gb, map_file, rm_mhc)


#########################################
### Significant interaction calling using Fit-Hi-C ###
#########################################
### create interaction count files
for chr in `seq 1 22`
do
cat RawCountMatrix40Kb/IMR90_40Kb_bin_chr"$chr"_raw_count_matrix.txt | awk -v chr=$chr '{start=(NR-1)*40000+20000; for (j=(NR+1);j<=(NR+75);j++){end=(j-1)*40000+20000; if (((end-start)<=3000000)&&(j<=NF)) print chr "\t" start "\t" chr "\t" end "\t" int($j+0) ; } }' >Fit_HiC_input_40Kb/IMR90.40Kb.raw.chr"$chr".txt
done

#### create marginal count files
for chr in `seq 1 22`
do
cat RawCountMatrix40Kb/IMR90_40Kb_bin_chr"$chr"_raw_count_matrix.txt | awk -v chr=$chr '{start=(NR-1)*40000+20000; suma=0; for (j=(NR-75);j<=(NR+75);j++){if ((j!=NR)&&(j>=1)&&(j<=NF)) suma=suma+$j;} print chr "\t0\t" start "\t" int(suma) "\t" 1; }' >Fit_HiC_input_40Kb/IMR90.40Kb.raw.marginal.chr"$chr".txt
done

gzip Fit_HiC_input_40Kb/*.txt

for chr in `seq 1 22`
do
fit-hi-c.py -L 20000 -U 3000000 -i Fit_HiC_input_40Kb/IMR90.40Kb.raw.chr"$chr".txt.gz -f Fit_HiC_input_40Kb/IMR90.40Kb.raw.marginal.chr"$chr".txt.gz -o Fit_HiC_output_40Kb/IMR90 -l fithic_40Kb_IMR90_"$chr"
done