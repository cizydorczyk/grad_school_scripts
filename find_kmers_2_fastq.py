from sys import argv
# Specify input fastqc data file, fastq file with reads, and output file
# in the command line
script, fastqcdatafile, fastqfile, outputfastqfile = argv

# function to get list of kmers from fastqc_data.txt file:
def get_kmers(fastqcdatafile):
    ff = open(fastqcdatafile, "r")
    list1 = list(ff)
    # for each indexed line in the file, if the line starts with "#Sequence", change the value of list1 to contain all
    # lines from that line's index+1 until the end of the file, minus the very last line
    kmer_list = []
    for lnum, line in enumerate(list1):
        if line.startswith("#Sequence"):
            list3 = list1[lnum+1:-1]
            for item in list3:
                list2 = item.split("\t")
                kmer_list.append(list2[0])
        else:
            list3 = []
    ff.close()
    return kmer_list

kmers = get_kmers(fastqcdatafile)
print kmers


# Create output file
outfile = open(outputfastqfile, 'w')
outfile.close()


# Open output file for appending
outfile = open(outputfastqfile, 'a')

# Open fastq file, read line by line, if kmer is in line, write line to
# output file, along with the header for that read, in fasta format

gg = open(fastqfile, 'r')
list2 = list(gg)
gg.close()

totalcount = 0
uniqcount = 0
list3 = []
for lnum, line in enumerate(list2):
    for i in kmers:
        if i in line:
            totalcount += 1
            if line not in list3:
                uniqcount += 1
                outfile.write(list2[lnum-1])
                outfile.write(list2[lnum])
                outfile.write(list2[lnum+1])
                outfile.write(list2[lnum+2])
                list3.append(line)
outfile.close()
print "Total reads with kmers: " + str(totalcount)
print "Unique reads with kmers: " + str(uniqcount)

print list3
