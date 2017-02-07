from sys import argv

script, inputfile = argv

with open(inputfile, 'r') as infile:
    list1 = list(infile)

dict1 = {}

for num, i in enumerate(list1):
    if i.startswith(">"):
        #print i,
        tick = 1
        sequence = ''
        j = list1[num+tick]
        while not j.startswith(">"):
            sequence += j.strip()
            tick += 1
            try:
                j = list1[num+tick]
            except IndexError:
                break
        dict1[i] = sequence


def GC_Content(sequence):
    count = 0.0
    for i in sequence:
        if i in "GCgc":
            count += 1.0
    return count/float(len(sequence))*100.00

for key, sequence in dict1.iteritems():
    print key,
    print GC_Content(sequence)
