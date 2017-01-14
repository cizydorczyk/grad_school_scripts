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

class FastqObject(object):
    def __init__(self, header, sequence, spacer, quality):
        self.header = header
        self.sequence = sequence
        self.spacer = spacer
        self.quality = quality

def Fastq_Parser(fastqfile):
    fastqobjects = []
    print "\tNow parsing: " + str(fastqfile)
    with open(fastqfile, 'r') as infile:
        reads = [line.strip("\n") for line in infile]
    count = 0
    for lnum, line in enumerate(reads):
        temp1 = []
        if line.startswith("@"):
            count += 1
            temp1.append(line.strip("\n"))
            temp1.append(reads[lnum+1])
            temp1.append(reads[lnum+2])
            temp1.append(reads[lnum+3])
            fastqobjects.append(FastqObject(temp1[0], temp1[1], temp1[2], temp1[3]))
    print "\tNumber of reads in file: " + str(count)
    print "\tNumber of fastq objects created: " + str(len(fastqobjects))
    return fastqobjects

a = Fastq_Parser(fastqfile)
print len(a)

def find_reads_with_kmers(fastqfileobjectlist, kmer_list):
    total_count = 0
    uniq_count = 0
    already_counted_reads = []
    with open(outputfastqfile, 'w') as outfile:
        for read in fastqfileobjectlist:
            for kmer in kmer_list:
                if kmer in read.sequence and read.header not in already_counted_reads:
                    total_count += 1
                    outfile.write(str(read.header) + "\n")
                    outfile.write(str(read.sequence) + "\n")
                    outfile.write(str(read.spacer) + "\n")
                    outfile.write(str(read.quality) + "\n")
                    already_counted_reads.append(read.header)
    print "Reads with kmers: " + str(total_count)

#def unique_reads(fastqfileobjectlist):
#    unique_reads = []
#    unique_count = 0
#    for read in fastqfileobjectlist:
#        if read.sequence not in unique_reads:
#            unique_count += 1
#            unique_reads.append(read.sequence)
#    print "Unique reads: " + str(unique_count)
#    print "Duplicate reads: " + str(len(fastqfileobjectlist) - unique_count)
b = find_reads_with_kmers(a, kmers)
#unique_reads(a)





# Create output file
#outfile = open(outputfastqfile, 'w')
#outfile.close()


# Open output file for appending
#outfile = open(outputfastqfile, 'a')

# Open fastq file, read line by line, if kmer is in line, write line to
# output file, along with the header for that read, in fasta format

#gg = open(fastqfile, 'r')
#list2 = list(gg)
#gg.close()

#totalcount = 0
#uniqcount = 0
#list3 = []
#for lnum, line in enumerate(list2):
#    for i in kmers:
#        if i in line:
#            totalcount += 1
#            if line not in list3:
#                uniqcount += 1
#                outfile.write(list2[lnum-1])
#                outfile.write(list2[lnum])
#                outfile.write(list2[lnum+1])
#                outfile.write(list2[lnum+2])
#                list3.append(line)
#outfile.close()
#print "Total reads with kmers: " + str(totalcount)
#print "Unique reads with kmers: " + str(uniqcount)
