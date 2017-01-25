from sys import argv


script, fastq1, fastq2, kmers1, kmers2 = argv

class FastqObject(object):
    def __init__(self, header, shortheader, sequence, spacer, quality, fullread):
        self.header = header
        self.shortheader = shortheader
        self.sequence = sequence
        self.spacer = spacer
        self.quality = quality
        self.fullread = fullread

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
            temp1.append(line)
            temp19 = line.split(' ')
            temp1.append(temp19[0])
            temp1.append(reads[lnum+1])
            temp1.append(reads[lnum+2])
            temp1.append(reads[lnum+3])
            temp1.append('\n'.join(temp1))
            fastqobjects.append(FastqObject(temp1[0], temp1[1], temp1[2], temp1[3], temp1[4], temp1[5]))
    print "\tNumber of reads in file: " + str(count)
    print "\tNumber of fastq objects created: " + str(len(fastqobjects))
    return fastqobjects

a = Fastq_Parser(fastq1)
b = Fastq_Parser(fastq2)
print "R1 reads total: " + str(len(a))
print "R2 reads total: " + str(len(b))
c = set(a)
d = set(b)
print "Length set a :" + str(len(c))
print "Length set b :" + str(len(d))

def reads_to_remove(kmerfile1, kmerfile2):
    with open(kmerfile1, 'r') as infile:
        headers1 = []
        for line in infile:
            if line.startswith("@"):
                temp1 = line.split(' ')
                headers1.append(temp1[0])
        #headers1 = [line.strip("\n") for line in infile if line.startswith("@")]
    with open(kmerfile2, 'r') as infile2:
        headers2 = [line.strip("\n") for line in infile2 if line.startswith("@")]
    headers = headers1 + headers2
    print len(headers1)
    print len(headers2)
    return headers
f = reads_to_remove(kmers1, kmers2)
g = set(f)
print len(g)


#def filter1(set_to_filter, headers_to_remove, outfastqfile):
#    count = 0
#    with open(outfastq, 'w') as outfile:
#        outfile.write('\n'.join(read.fullread for read in set_to_filter if read.header not in headers_to_remove))

#d = filter1(c, e, outfastq)
