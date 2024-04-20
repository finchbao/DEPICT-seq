#!/bin/bash
map_dir=~/DEPICT_seq/mapped/
cd ${map_dir}
output_dir=~/DEPICT_seq/qualimap_output/
for i in `ls *Aligned.sortedByCoord.out.bam`
do 
	i=${i/Aligned.sortedByCoord.out.bam} 
	qualimap rnaseq -bam ${map_dir}${i}Aligned.sortedByCoord.out.bam \
   -outdir ${output_dir}${i}/ \
   -gtf ~/DEPICT_seq/dataset/gencode.v41.annotation.gtf \
   -outformat PDF:HTML \
   --java-mem-size=10G \
	done
