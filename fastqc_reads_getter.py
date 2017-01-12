from sys import argv
import os.path

script, fastqc_data_file, outputtsv = argv

def get_reads(fastqcdatafile):
    with open(fastqc_data_file, 'r') as infile:
        for line in infile:
            if line.startswith("Filename"):
                a = line.strip("\n").split("\t")
            elif line.startswith("Total Sequences"):
                b = line.strip("\n").split("\t")
        total = a + b
    if os.path.isfile(outputtsv):
        with open(outputtsv, 'a') as outfile:
            outfile.write(total[1]+"\t")
            outfile.write(total[3]+"\n")
    elif not os.path.isfile(outputtsv):
        with open(outputtsv, 'w') as outfile:
            outfile.write(total[1]+"\t")
            outfile.write(total[3]+"\n")

get_reads(fastqc_data_file)
