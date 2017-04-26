#! /bin/bash

## Input: 
##	isolate number,	1
##	input_bwa_sam_file_directory,	2
##	output_bwa_sorted_indexed_bam_directory,	3
##	input_novo_sam_directory,	4
##	output_novo_sorted_indx_bam_directory,	5
##	input_last_sam_directory,	6
##	output_last_sorted_indexed_bam_directory,	7
##	bwa samtools indexed reference	8
##	novoalign samtools indexed reference	9
##	lastalign samtools indexed reference	10
##	bwa_vcf_output ddirectory	11
##	novo vcf output directory	12
##	last vcf output directory	13

~/Software/samtools-1.3.1/samtools view -bhuS $2$1"_aln.sam" | ~/Software/samtools-1.3.1/samtools sort -@ 8 -o $3$1"_aln.bam"
~/Software/samtools-1.3.1/samtools index $3$1"_aln.bam"
~/Software/samtools-1.3.1/samtools mpileup -ugf $8 $3$1"_aln.bam" | ~/Software/bcftools-1.3.1/bcftools call -m -v --ploidy 1 --threads 8 -o ${11}$1"_bwa.vcf"

~/Software/samtools-1.3.1/samtools view -bhuS $4$1"_aln.sam" | ~/Software/samtools-1.3.1/samtools sort -@ 8 -o $5$1"_aln.bam"
~/Software/samtools-1.3.1/samtools index $5$1"_aln.bam"
~/Software/samtools-1.3.1/samtools mpileup -ugf $9 $5$1"_aln.bam" | ~/Software/bcftools-1.3.1/bcftools call -m -v --ploidy 1 --threads 8 -o ${12}$1"_novo.vcf"

~/Software/samtools-1.3.1/samtools view -bhuS $6$1"_aln.sam" | ~/Software/samtools-1.3.1/samtools sort -@ 8 -o $7$1"_aln.bam"
~/Software/samtools-1.3.1/samtools index $7$1"_aln.bam"
~/Software/samtools-1.3.1/samtools mpileup -ugf ${10} $7$1"_aln.bam" | ~/Software/bcftools-1.3.1/bcftools call -m -v --ploidy 1 --threads 8 -o ${13}$1"_last.vcf"


