from sys import argv

script, file1, file2, out1, out2 = argv

# Read object class for each read in a fastq file
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

a = Fastq_Parser(file1)
b = Fastq_Parser(file2)

with open(out1, 'w') as outfile1:
    for i in a:
        outfile1.write(i.shortheader + '\n')
        outfile1.write(i.sequence + '\n')
        outfile1.write(i.spacer + '\n')
        outfile1.write(i.quality + '\n')

with open(out2, 'w') as outfile2:
    for i in b:
        outfile2.write(i.shortheader + '\n')
        outfile2.write(i.sequence + '\n')
        outfile2.write(i.spacer + '\n')
        outfile2.write(i.quality + '\n')
