from sys import argv

script, isolate_list, outputdir = argv

isolatelist = []
with open(isolate_list, 'r') as infile1:
    for line in infile1:
        isolatelist.append(line.strip())

for i in isolatelist:
    with open(outputdir + "/" + i + ".job", 'w+') as outfile1:
        outfile1.write('#!/bin/bash' +'\n' + '#PBS -l nodes=1:ppn=8' + '\n' +
        '#PBS -l walltime=48:00:00' + '\n' + '#PBS -N ' + i + '_pae' + '\n' +
        '#PBS -A zke-503-ab' + '\n' + '\n' + 'cd /scratch/d/dguttman/cizydor/'
        + '\n' + '\n' + 'module load blast/2.2.30+' + '\n' + '\n' + 'blastp -query /scratch/d/dguttman/cizydor/3.blastquery/' +
        i + '_pae.fasta' +
        ' -db /scratch/d/dguttman/cizydor/3.blastdb/goodProteins.fasta -seg yes -dbsize 100000000 -evalue 1e-05 -outfmt 6 -num_threads 8 -out /scratch/d/dguttman/cizydor/3.blastres/' +
        i + '_pae.tab')
