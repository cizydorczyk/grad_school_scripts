from sys import argv

script, outfile = argv


with open(outfile, 'w') as outfile1:

    outfile1.write("#!/bin/bash" + '\n')
    outfile1.write("#PBS -l nodes=1:ppn=8" + '\n')
    outfile1.write("#PBS -l walltime=00:30:00" + '\n')
    outfile1.write("#PBS -N " + str(outfile) + '\n')
    outfile1.write("#PBS -A zke-503-ab" + '\n' + '\n')
    outfile1.write("cd /scratch/d/dguttman/cizydor/" + '\n' + '\n')
    outfile1.write('module load blast/2.2.30+' + '\n' + '\n')
    outfile1.write("blastn -query /scratch/d/dguttman/cizydor/E4_assembly_contigs/" + str(outfile) + "_assembly.fa -db /scratch/d/dguttman/cizydor/blastdb_nt/blastdb/nt -evalue 1e-5 -out /scratch/d/dguttman/cizydor/blast_output/" + str(outfile) + "_blast_output -outfmt '7 qseqid sseqid pident length evalue sscinames' -num_threads 8")
