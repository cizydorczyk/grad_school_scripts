from sys import argv

script, fastqfile = argv

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

"""Returns the proper number of lines (same as grep) that start with @ in the fastq file"""
