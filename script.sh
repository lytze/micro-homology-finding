## 1. Making alignment
##    Because the following MTD read catching methods depends on clipped sequence, here we sed -Y on for bwa mem
bwa index ref/lib1.fa
bwa mem -t 8 -Y ref/lib1.fa data/Library_1_clean_R1.fastq.gz data/Library_1_clean_R2.fastq.gz | 
	samtools view -@ 2 -F 2048 -Sb |
	samtools sort -@ 2 -m 10G -o mapped/lib1.sort.bam

## 2. MH pair finding and generate signature features for MTD reads
code/find_mh probes/lib1.mh < ref/lib1.fa
code/generate_signatures.awk ref/lib1.fa \
    probes/lib1.mh.1.seq1.out probes/lib1.mh.2.seq2.out > probes/signatures.lib1.out

## 3. Count for read depth for MH positions
##    Because we used super-deep sequencing, this step took a long time
##    Run step 4 meanwhile to save time
samtools index -@ 4 mapped/lib1.sort.bam
code/sigtobed.awk probes/signatures.lib1.out > probes/lib1.pos.bed
samtools depth -b probes/lib1.pos.bed -d 0 mapped/lib1.sort.bam > mapped/lib1.cov

## 4. Catch Reads with MTD signature features from the alignment map
samtools view mapped/lib1.sort.bam | code/catch_signatures.awk probes/signatures.lib1.out > results/lib1.sign.count.tsv

## 5. Normalize the count number based on read-depths
code/normalize_count.awk mapped/lib1.cov results/lib1.sign.count.tsv