from sys import argv
import os.path

script, fastqfile, outputtsv = argv

def get_total_bases_from_fastq(fastqfile):
    print "Parsing "+fastqfile
    bases_count = 0
    with open(fastqfile, 'r') as infile:
        for line in infile:
            if line.startswith("@"):
                bases_count += len(next(infile).strip())
    if os.path.isfile(outputtsv):
        with open(outputtsv, 'a') as outfile:
            outfile.write(fastqfile + "\t")
            outfile.write(str(bases_count)+"\n")
    elif not os.path.isfile(outputtsv):
        with open(outputtsv, 'w') as outfile:
            outfile.write(fastqfile + "\t")
            outfile.write(str(bases_count)+"\n")


get_total_bases_from_fastq(fastqfile)
