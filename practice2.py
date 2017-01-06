from sys import argv
script, fastqcfile, fastqfile, outputfastafile = argv

# function to get list of kmers from fastqc_data.txt file:
def get_kmers(fastqcdatafile):
    ff = open(fastqcdatafile, "r")
    list1 = list(ff)
    # for each indexed line in the file, if the line starts with "#Sequence", change the value of list1 to contain all
    # lines from that line's index+1 until the end of the file, minus the very last line
    for lnum, line in enumerate(list1):
        if line.startswith("#Sequence"):
            list1 = list1[lnum+1:-1]
    ff.close()

    kmer_list = []
    for item in list1:
        list2 = item.split("\t")
        kmer_list.append(list2[0])
    return kmer_list

kmers = get_kmers(fastqcfile)
print kmers

outfile = open(outputfastafile, 'w')
outfile.close()

outfile = open(outputfastafile, 'a')
prevline = ""

with open(fastqfile, 'r') as infile:
    for line in infile:
        for i in kmers:
            if i in line:
                outfile.write(">"+prevline)
                outfile.write(line)
outfile.close()
