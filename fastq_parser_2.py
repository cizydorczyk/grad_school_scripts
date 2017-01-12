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
    print "Now parsing: " + str(fastqfile)
    with open(fastqfile, 'r') as infile:
        numline = 0
        for line in infile:
            if line.startswith("@"):
                numline += 1
    return numline

print Fastq_Parser(fastqfile)

"""Returns the proper number of lines (same as grep) that start with @ in the fastq file"""
