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

fastq_file_reads = Fastq_Parser(fastq_file_to_filter)

def reads_to_remove(toremove_fastq):
    print "\tIdentifying reads to remove in: " + str(toremovefastq)
    readstoremove = []
    with open(toremove_fastq, 'r') as infile:
        for line in infile:
            if line.startswith("@"):
                line = next(infile)
                readstoremove.append(line.strip("\n"))
    print "\tNumber of reads in kmer file: " + str(len(readstoremove))
    return readstoremove

readstoremove = reads_to_remove(toremovefastq)


def remove_reads(fastq_reads_objects_list, list_of_reads_to_remove):
    print "Removing reads from: " + fastq_file_to_filter
    count = 0
    for i in fastq_reads_objects_list:
        if i.sequence in list_of_reads_to_remove:
            count += 1
    return count

temp1 = remove_reads(fastq_file_reads, readstoremove)
print temp1
