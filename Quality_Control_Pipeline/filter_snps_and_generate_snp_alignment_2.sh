#!/bin/bash

########### REMEMBER TO SET REF GENOME LENGTH - 150 bp BELOW!!! #################

# Set reference genome name:
REFGENOME=PAO1

# Set tool paths (generally don't change):

BCFTOOLS=/home/conrad/Software/bcftools-1.3.1/bcftools
BGZIP=/home/conrad/Software/htslib-1.5/bgzip
TABIX=/home/conrad/Software/htslib-1.5/tabix
SNPSITES=/usr/bin/snp-sites

# Set python script paths (generally don't change):
CALLSNPS=/home/conrad/grad_school_scripts/Quality_Control_Pipeline/call_snps_on_vcf_1.py
COMBINEFASTASEQ=/home/conrad/grad_school_scripts/Quality_Control_Pipeline/combine_fasta_sequences_1.py
GETISOLATESNPS=/home/conrad/grad_school_scripts/Quality_Control_Pipeline/get_isolate_snps_1.py
ANNOTATESNPS=/home/conrad/grad_school_scripts/Quality_Control_Pipeline/annotation_script_4.py

# Set project folder:

PF=/home/conrad/Data/primary_project_3/reference_alignments/G49_PAO1

# Set printed-to-screen script text color (generally don't change):
COLOR='\033[1;36m'
NC='\033[0m'

# Input variables:

# 1 = reference genome (fasta file)
# 2 = isolate list (text file, one isolate number per line)
# 3 = TAB genome annotation (from pseudomonas.com) - make sure there are no blank fields in the file!!
# 4 = GTF genome annotation (from pseudomonas.com) - make sure there are no blank fields in the file!!

cd $PF

mkdir uncalled_filtered_vcf_files

cd raw_vcf_files

printf "${COLOR}Filtering raw vcf files...${NC}\n"
################ SET REF LEN - 150 bp IN LINE BELOW WHERE POS<(some number) #######################
for i in $(cat $2); do $BCFTOOLS filter -g 150 -i 'QUAL>=50 & DP>=20 & TYPE="snp" & POS>150 & POS<6264254 & DP4[2]>3 & DP4[3]>3 & MQ>25' -o $PF/uncalled_filtered_vcf_files/$i"_uncalled_filtered_variants.vcf" $i"_raw_variants.vcf"; done

cd ..

mkdir snps_called_filtered_vcf_files

cd snps_called_filtered_vcf_files

printf "${COLOR}Calling snps on filtered vcf files...${NC}\n"
for i in $(cat $2); do python $CALLSNPS $PF/uncalled_filtered_vcf_files/$i"_uncalled_filtered_variants.vcf" $i"_snps_called_filtered_snps.vcf"; done

cd ..

mkdir whole_genome_consensus_sequences

cp $1 $PF/whole_genome_consensus_sequences/reference_wg.fa

cd snps_called_filtered_vcf_files

printf "${COLOR}Generating consensus whole genome sequences...${NC}\n"
for i in $(cat $2); do $BGZIP $i"_snps_called_filtered_snps.vcf"; $TABIX $i"_snps_called_filtered_snps.vcf.gz"; $BCFTOOLS consensus -f $1 $i"_snps_called_filtered_snps.vcf.gz" > $PF/whole_genome_consensus_sequences/$i"_wg_consensus.fa"; done

cd ..

mkdir alignments

cd whole_genome_consensus_sequences

printf "${COLOR}Generating consensus whole genome alignment...${NC}\n"
for i in *.fa; do python $COMBINEFASTASEQ $i $PF/alignments/whole_genome_consensus_alignment.fa; done

cd ../alignments

printf "${COLOR}Generating SNP alignment...${NC}\n"
$SNPSITES -o snp_alignment.fa $PF/alignments/whole_genome_consensus_alignment.fa

cd ../snps_called_filtered_vcf_files

for i in *.gz; do $BGZIP -d $i; done

for i in *.vcf; do realpath $i >> snps_called_filtered_vcf_file_list.txt; done

cd ..

mkdir unannotated_snps

printf "${COLOR}Generating list of SNPs for each isolate...${NC}\n"
python $GETISOLATESNPS $PF/snps_called_filtered_vcf_files/snps_called_filtered_vcf_file_list.txt $PF/unannotated_snps/

mkdir annotated_snps

cp $3 $PF/reference/$REFGENOME"_TAB_annotation.txt"
cp $4 $PF/reference/$REFGENOME"_GTF_annotation.txt"

cd unannotated_snps

printf "${COLOR}Annotating SNPs for each isolate...${NC}\n"
for i in $(cat $2); do python $ANNOTATESNPS $3 $4 $i"_unannotated_snps.txt" $PF/annotated_snps/$i"_annotated_snps.txt"; done
