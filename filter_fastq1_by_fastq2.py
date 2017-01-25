from sys import argv


script, fastq1, fastq2, outfastq = argv

class FastqObject(object):
    def __init__(self, header, sequence, spacer, quality, fullread):
        self.header = header
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
            temp1.append(reads[lnum+1])
            temp1.append(reads[lnum+2])
            temp1.append(reads[lnum+3])
            temp1.append('\n'.join(temp1))
            fastqobjects.append(FastqObject(temp1[0], temp1[1], temp1[2], temp1[3], temp1[4]))
    print "\tNumber of reads in file: " + str(count)
    print "\tNumber of fastq objects created: " + str(len(fastqobjects))
    return fastqobjects

a = Fastq_Parser(fastq2)
c = set(a)

def reads_to_remove(fastqfile):
    print "Identifying reads in " + fastq1 + " to remove"
    with open(fastqfile, 'r') as infile:
        reads = [line.strip("\n") for line in infile if line.startswith("@")]
    print "Found " + str(len(reads)) + " reads to remove"
    return reads
b = reads_to_remove(fastq1)
e = set(b)


def filter1(set_to_filter, headers_to_remove, outfastqfile):
    count = 0
    with open(outfastq, 'w') as outfile:
        outfile.write('\n'.join(read.fullread for read in set_to_filter if read.header not in headers_to_remove))

#d = filter1(c, e, outfastq)
