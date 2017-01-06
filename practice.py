f = "/Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_output_unzipped/6_R1_fastqc/fastqc_data.txt"

# function to get list of kmers from fastqc_data.txt file:
def get_kmers(fastqcdatafile):
    ff = open(fastqcdatafile, "r")
    list1 = list(ff)
    # for each indexed line in the file, if the line starts with "#Sequence", change the value of list1 to contain all
    # lines from that line's index+1 until the end of the file, minus the very last line
    for lnum, line in enumerate(list1):
        if line.startswith("#Sequence"):
            list1 = list1[lnum+1:-1]
    ff.close()

    kmer_list = []
    for item in list1:
        list2 = item.split("\t")
        kmer_list.append(list2[0])
    return kmer_list

kmers = get_kmers(f)

g = "/Users/conradizydorczyk/Dropbox/6_R1_trimmed.fastq"

h = open("/Users/conradizydorczyk/Desktop/testfile2.txt", 'w')
h.close()

h = open("/Users/conradizydorczyk/Desktop/testfile2.txt", 'a')

prevline = ""

for i in kmers:
    with open(g, 'r') as fg:
        for line in fg:
            if i in line:
                h.write(prevline)
                h.write(line)
            prevline = line
h.close()