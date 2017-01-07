from sys import argv
# Specify input fastqc data file, fastq file with reads, and output file
# in the command line
script, fastqcdatafile, fastqfile, outputfastafile = argv

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

kmers = get_kmers(fastqcdatafile)
print kmers

# Create output file
outfile = open(outputfastafile, 'w')
outfile.close()


# Open output file for appending
outfile = open(outputfastafile, 'a')

# Open fastq file, read line by line, if kmer is in line, write line to
# output file, along with the header for that read, in fasta format

gg = open(fastqfile, 'r')
list2 = list(gg)
gg.close()

count = 0
list3 = []
for lnum, line in enumerate(list2):
    for i in kmers:
        if i in line and line not in list3:
            count += 1
            outfile.write(list2[lnum-1])
            outfile.write(list2[lnum])
            outfile.write(list2[lnum+1])
            outfile.write(list2[lnum+2])
            list3.append(line)
outfile.close()
print "Reads with kmers:", count
