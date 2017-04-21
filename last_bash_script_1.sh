#! /bin/sh

# last_bash_script_1.sh <fastq file directory> <output directory> <isolate number> <full path to reference prefix>

~/Software/last-849/scripts/fastq-interleave $1$3_R1.fastq $1/$3_R2.fastq > $2$3_interleaved.fastq

~/Software/last-849/src/lastal -v -P0 -Q1 -i1 $4 $2$3_interleaved.fastq > $2$3_temp.maf

~/Software/last-849/src/last-pair-probs $2$3_temp.maf > $2$3_final.maf

~/Software/last-849/scripts/maf-convert sam $2$3_final.maf > $2$3_aln.sam

rm $2$3_temp.maf
rm $2$3_interleaved.fastq
rm $2$3_final.maf
