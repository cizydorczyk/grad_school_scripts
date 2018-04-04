#!/bin/bash

########### REMEMBER TO SET REF GENOME LENGTH - 150 bp BELOW!!! #################

# Set reference genome name:
REFGENOME=LESB58

# Set tool paths (generally don't change):
SAMTOOLS=/home/conrad/Software/samtools-1.3.1/samtools
PICARD=/home/conrad/Software/picard.jar
GATK=/home/conrad/Software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
BCFTOOLS=/home/conrad/Software/bcftools-1.3.1/bcftools

# Set python script paths (generally don't change):
INDELLIST=/home/conrad/grad_school_scripts/Quality_Control_Pipeline/indel_list_1.py
ANNOTATEINDELS=/home/conrad/grad_school_scripts/Quality_Control_Pipeline/indel_annotation_script_2.py

# Set project folder:
PF=/home/conrad/Data/primary_project_3/reference_alignments/H18_lesb58_indels

# Set directory with sorted, indexed bam files with read groups added (must do this prior to running this script):
BAMINPUT=/home/conrad/Data/primary_project_3/reference_alignments/H18_lesb58_indels/sorted_indexed_bam_files_with_rg

# Set printed-to-screen script text color (generally don't change):
COLOR='\033[1;36m'
NC='\033[0m'

# Input variables:

# 1 = reference genome (fasta file)
# 2 = isolate list (text file, one isolate number per line)
# 3 = TAB genome annotation (from pseudomonas.com) - make sure there are no blank fields in the file!!
# 4 = GTF genome annotation (from pseudomonas.com) - make sure there are no blank fields in the file!!

cd $PF

# mkdir reference_indels

# cp $1 $PF/reference_indels/$REFGENOME".fasta"
# cp $3 $PF/reference_indels/$REFGENOME"_TAB_annotation_no_blanks.txt"
# cp $4 $PF/reference_indels/$REFGENOME"_GTF_annotation_no_blanks.txt"

cd reference

# $SAMTOOLS faidx $REFGENOME".fasta"
# java -jar $PICARD CreateSequenceDictionary REFERENCE=$REFGENOME".fasta" OUTPUT=$REFGENOME".dict"
java -jar $PICARD CreateSequenceDictionary REFERENCE=$1 OUTPUT="reference.dict"
cd ..

mkdir raw_vcf_files

cd raw_vcf_files

mkdir picard_markduplicate_metrics

for i in $(cat $2); do java -Xmx20g -jar $PICARD MarkDuplicates I=$BAMINPUT/$i"_aln_rg.bam" O=$i"_aln_rg_markeddup.bam" M=picard_markduplicate_metrics/$i"_marked_duplicate_metrics.txt"; done

for i in $(cat $2); do $SAMTOOLS index $i"_aln_rg_markeddup.bam"; done

for i in $(cat $2); do java -Xmx20g -jar $GATK -T HaplotypeCaller -R $PF/reference/$REFGENOME".fasta" -I $i"_aln_rg_markeddup.bam" -o $i"_raw.vcf" -ploidy 1 -nct 8 -stand_call_conf 20; done

cd ..

mkdir filtered_indel_vcf_files

cd filtered_indel_vcf_files

for i in $(cat $2); do java -Xmx20g -jar $GATK -T SelectVariants -R $PF/reference/$REFGENOME".fasta" -nt 8 -V $PF/raw_vcf_files/$i"_raw.vcf" -selectType INDEL -o $i"_raw_indels.vcf"; done

################ SET REF LEN - 150 bp IN LINE BELOW WHERE POS<(some number) #######################
for i in $(cat $2); do $BCFTOOLS filter -i 'QUAL>=30 & DP>=20 & POS>150 & POS<6264254' -o $i"_filtered_indels.vcf" --threads 8 $PF/filtered_indel_vcf_files/$i"_raw_indels.vcf"; done

cd ..

mkdir unannotated_indels
mkdir annotated_indels

for i in $(cat $2); do python $INDELLIST $PF/filtered_indel_vcf_files/$i"_filtered_indels.vcf" $PF/unannotated_indels/$i"_unannotated_indels.txt"; done
for i in  $(cat $2); do python $ANNOTATEINDELS $3 $4 $PF/unannotated_indels/$i"_unannotated_indels.txt" $PF/annotated_indels/$i"_annotated_indels.txt"; done

cd raw_vcf_files

for i in $(cat $2); do rm $i"_aln_rg_markeddup.bam"; rm $i"_aln_rg_markeddup.bam.bai"; rm $i"_raw.vcf.idx"; done

cd ../filtered_indel_vcf_files

for i in $(cat $2); do rm $i"_raw_indels.vcf"; rm $i"_raw_indels.vcf.idx"; done
