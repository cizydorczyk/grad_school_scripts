from sys import argv
import os.path

script, fastqc_data_file, outputtsv = argv

def get_avg_length(fastqcdatafile):
    print "Parsing "+fastqc_data_file
    with open(fastqc_data_file, 'r') as infile:
        for line in infile:
            if line.startswith("Filename"):
                a = line.strip("\n").split("\t")
                filename = a[1]
            elif line.startswith("#Length"):
                totalbases = 0
                line = next(infile)
                while not line.startswith(">>"):
                    b = line.strip("\n").split("\t")
                    totalbases += float(b[0])*float(b[1])
                    line = next(infile)
    if os.path.isfile(outputtsv):
        with open(outputtsv, 'a') as outfile:
            outfile.write(filename + "\t")
            outfile.write(str(totalbases)+"\n")
    elif not os.path.isfile(outputtsv):
        with open(outputtsv, 'w') as outfile:
            outfile.write(filename + "\t")
            outfile.write(str(totalbases)+"\n")
get_avg_length(fastqc_data_file)
