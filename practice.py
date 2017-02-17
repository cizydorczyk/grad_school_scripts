from sys import argv
from oset import oset
script, fastq1, contaminants = argv

class FastqObject(object):
    def __init__(self, shortheader, header, sequence, spacer, quality, fullread):
        self.shortheader = shortheader
        self.header = header
        self.sequence = sequence
        self.spacer = spacer
        self.quality = quality
        self.fullread = fullread

# Fastq_Parser parses a fastq file and calls the FastqObject class on each read,
# creating a unique read object for each read in the file
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
            temp19 = line.split(' ')
            temp1.append(temp19[0])
            temp1.append(line)
            temp1.append(reads[lnum+1])
            temp1.append(reads[lnum+2])
            temp1.append(reads[lnum+3])
            temp1.append('\n'.join(temp1[1:]))
            fastqobjects.append(FastqObject(temp1[0], temp1[1], temp1[2], temp1[3], temp1[4], temp1[5]))
    print "\tNumber of reads in file: " + str(count)
    print "\tNumber of fastq objects created: " + str(len(fastqobjects))
    return fastqobjects

a = Fastq_Parser(fastq1)
#b = Fastq_Parser(fastq2)
print "R1 reads total: " + str(len(a))
#print "R2 reads total: " + str(len(b))
# Convert list outputs from Fastq_Parser for each of the two input fastq files
# into ordered sets, which are much faster to test for membership in than lists (used
# later)
c = oset(a)
#d = oset(b)
print "Length set a :" + str(len(c))
#print "Length set b :" + str(len(d))

# Creates a list of all kmer headers (minus the #:#:#:# bit at the end of each header)
# from both input kmer files for both input (R1 and R2) fastq files. This list
# may contain duplicates at this point, but they are removed when the output
# list is converted into a set (see below)
def reads_to_remove(kmerfile1):
    with open(contaminants, 'r') as infile:
        headers1 = []
        for line in infile:
            if line.startswith("@"):
                temp1 = line.split(' ')
                headers1.append(temp1[0])
        #headers1 = [line.strip("\n") for line in infile if line.startswith("@")]
    headers = headers1
    print "Number of kmers in file1: " + str(len(headers1))
    #print "Number of kmers in file2: " + str(len(headers2))
    return headers
f = reads_to_remove(contaminants)

# Convert output from reads_to_remove into an ordered set.
# Sets only contain UNIQUE entries; so if two short headers are identical, they
# will only appear once in the set (hence the length of set g may be lower
# than the sum of the length of the headers list)
g = oset(f)
print "Number of unique kmer headers: " + str(len(g))
