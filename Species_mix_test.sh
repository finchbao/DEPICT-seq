#!/bin/bash
# Map the reads into hsa mm mix genome ref
trim_dir=~/DEPICT_seq/trim/
hm_map_dir=~/DEPICT_seq/hm_mapped/
hm_map_counts_dir=~/DEPICT_seq/hm_mapped_output/
cd ${trim_dir}
for file in $(ls -F |grep R1_001_val_1.fq.gz)
do 
	file1=${file/R1_001.fastq.gz}
  STAR --runThreadN 10  \
     --runMode alignReads  \
     --readFilesCommand zcat \
     --twopassMode Basic \
     --outSAMtype BAM SortedByCoordinate \
     --genomeDir ~/DEPICT_seq/GRCh38_and_GRCm39_index/  \
     --readFilesIn ${trim_dir}${file1}_R1_001_val_1.fq.gz ${trim_dir}${file1}_R2_001_val_2.fq.gz  \
     --outFileNamePrefix ${hm_map_dir}${file1}
	done
# counts the reads number
cd ${hm_map_dir}
for i in `ls *Aligned.sortedByCoord.out.bam`
do 
	i=${i/Aligned.sortedByCoord.out.bam}
	samtools index $hm_map_dir}${i}Aligned.sortedByCoord.out.bam
	samtools idxstats ${hm_map_dir}${i}Aligned.sortedByCoord.out.bam | awk '{print $1" "$3}' > ${hm_map_counts_dir}${i}.txt
	done
