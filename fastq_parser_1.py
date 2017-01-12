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
                features = []
                features.append(line.strip("\n"))
                line = next(infile)
                while not line.startswith("@"):
                    features.append(line.strip("\n"))
                    line = next(infile)
                fastqobjects.append(FastqObject(features[0], features[1], features[2], features[3]))
    print "Number of sequences in fastq file: " + str(numline)
    print "Number of objects in fastq file: " + str(len(fastqobjects))
    return fastqobjects
list_of_fastq_objects = Fastq_Parser(fastqfile)
