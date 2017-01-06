#!/bin/bash

for i in $(cat /Users/conradizydorczyk/grad_school_data/primary_project/file_list.txt)

do
	/Users/conradizydorczyk/grad_school_data/scripts/FastqcSectionGetter.py -s quality -i /Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_output_unzipped/$i"_R1_fastqc"/fastqc_data.txt -o /Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqcsectiongetter_quality/$i"_R1_quality.txt"
	/Users/conradizydorczyk/grad_school_data/scripts/FastqcSectionGetter.py -s quality -i /Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_output_unzipped/$i"_R2_fastqc"/fastqc_data.txt -o /Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqcsectiongetter_quality/$i"_R2_quality.txt"
done


#/home/conrad/My_Scripts/FastqcSectionGetter.py -s quality -i fastqc_data.txt -v works to run the command...
#Outputting a directory doesn't work...
#inputting a path works...
#Specifying path to specicif output file works...
#Essentially loop twice...to specify names properly...I ran it in 2 separate loops, but Julio told me to just put the 2nd run into the first loop...hopefully it works but I never tried.
