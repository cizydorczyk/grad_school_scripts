#!/bin/bash

for i in $(cat /home/conrad/Data/160902_NextSeq/primary_project/file_list.txt)

do
	/home/conrad/Software/FastQC/fastqc -o /home/conrad/Data/160902_NextSeq/primary_project/Fastqc/Fastqc_raw $i"_R1.fastq" $i"_R2.fastq"

done
