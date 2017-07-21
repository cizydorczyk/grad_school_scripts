from sys import argv
import os.path

script, fastqfile, outputtsv = argv

def get_total_bases_from_fastq(fastqfile):
    print "Parsing "+fastqfile
    read_count = 0
    with open(fastqfile, 'r') as infile:
        for line in infile:
            if line.startswith("@"):
                read_count += 1
    if os.path.isfile(outputtsv):
        with open(outputtsv, 'a') as outfile:
            outfile.write(fastqfile + "\t")
            outfile.write(str(read_count)+"\n")
    elif not os.path.isfile(outputtsv):
        with open(outputtsv, 'w') as outfile:
            outfile.write(fastqfile + "\t")
            outfile.write(str(read_count)+"\n")


get_total_bases_from_fastq(fastqfile)
