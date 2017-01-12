from sys import argv

script, toremovefastq, fastq_file_to_filter = argv

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
    print "\tNumber of sequences in fastq file: " + str(numline)
    print "\tNumber of objects in fastq file: " + str(len(fastqobjects))
    return fastqobjects

fastq_file_reads = Fastq_Parser(fastq_file_to_filter)

def reads_to_remove(toremove_fastq):
    print "\tIdentifying reads to remove in: " + str(toremovefastq)
    readstoremove = []
    with open(toremove_fastq, 'r') as infile:
        for line in infile:
            if line.startswith("@"):
                line = next(infile)
                readstoremove.append(line.strip("\n"))
    print "\tNumber of reads to remove identified: " + str(len(readstoremove))
    return readstoremove

readstoremove = reads_to_remove(toremovefastq)


def remove_reads(fastq_reads_objects_list, list_of_reads_to_remove):
    count = 0
    for i in fastq_reads_objects_list:
        if i.sequence in list_of_reads_to_remove:
            count += 1
    return count

temp1 = remove_reads(fastq_file_reads, readstoremove)
print temp1
