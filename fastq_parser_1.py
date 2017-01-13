from sys import argv

script, fastqfile = argv

class FastqObject(object):
    def __init__(self, header):
        self.header = header
        #self.sequence = sequence
        #self.spacer = spacer
        #self.quality = quality

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
#                nextline = next(infile)
#                while not nextline.startswith("@"):
#                    features.append(nextline.strip("\n"))
#                    nextline = next(infile)
                fastqobjects.append(FastqObject(features[0]))
    print "Number of sequences in fastq file: " + str(numline)
    print "Number of objects in fastq file: " + str(len(fastqobjects))
    return fastqobjects
list_of_fastq_objects = Fastq_Parser(fastqfile)

"""Returns same # of lines as grep version 2 of this script that start with @ when I comment out the above sections..."""
