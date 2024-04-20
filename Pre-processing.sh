#!/bin/bash
# Trim and mapping with STAR
workdir=~/DEPICT_seq/
trim_dir=~/DEPICT_seq/trim/
map_dir=~/DEPICT_seq/mapped/
for file in $(ls -F |grep R1_001.fastq.gz)
do 
	file1=${file/R1_001.fastq.gz}
	trim_galore -q 20 \
     --phred33 \
     --stringency 3 \
     --length 20 \
     -e 0.1 \
     --paired ${workdir}/${file1}R1_001.fastq.gz ${workdir}/${file1}R2_001.fastq.gz \
     --gzip -o ${trim_dir}
  STAR --runThreadN 10  \
     --runMode alignReads  \
     --readFilesCommand zcat \
     --twopassMode Basic \
     --outSAMtype BAM SortedByCoordinate \
     --genomeDir ~/DEPICT_seq/hg38_index/  \
     --readFilesIn ${trim_dir}${file1}_R1_001_val_1.fq.gz ${trim_dir}${file1}_R2_001_val_2.fq.gz  \
     --outFileNamePrefix ${map_dir}${file1}
	done

 # Counts the expression matrix
featureCounts -t gene \
  -g gene_id  \
  -a /data/baokaixuan/JZ_lab/20220831Flash-seq/dataset/gencode.v41.annotation.gtf \
  -o ${workdir}depict_mtx.txt  ${map_dir}*.sortedByCoord.out.bam
