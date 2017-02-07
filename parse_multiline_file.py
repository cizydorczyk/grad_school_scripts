# Eg. for multiline fasta or fastq...also includes method for when end of file
# is reached (except IndexError part)

from sys import argv

script, inputfile = argv

with open(inputfile, 'r') as infile:
    list1 = list(infile)

for num, i in enumerate(list1):
    if i.startswith(">"):
        print i,
        tick = 1
        sequence = ''
        i = list1[num+tick]
        while not i.startswith(">"):
            sequence += i.strip()
            tick += 1
            try:
                i = list1[num+tick]
            except IndexError:
                break
        print sequence