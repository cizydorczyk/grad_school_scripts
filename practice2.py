from sys import argv
import re

script, inputfile = argv

dict1 = {}
#should_print = False
with open(inputfile, 'r') as infile:
    file_list = list(infile)
for lnum, line in enumerate(file_list):
    if line.startswith(">"):
        print line
        i = 1
        while ">" not in line:
            for lnum, line in enumerate(file_list):
                print file_list[lnum+i]
                i += 1
                line = file_list[lnum+i]




            # should_print becomes True if was False and becomes False if was True
        #    should_print = not should_print
        #if should_print:
        #    print(line)
            # keep printing until next ">" encountered; stop at that point
            # and move on to next line that starts with ">"...



#        if line.startswith(">"):
#            header = line
#            print header
#print dict1

            #dict1[line.strip("\n")] = next(infile).strip("\n")

#def GC(DNA):
#    count = 0
#    for i in DNA:
#    print count
#    print len(DNA)
#    count = float(count)/float(len(DNA)) * 100.00
#    return count
#for key, value in dict1.iteritems():
#    print key
#    print GC(value)

#print len("CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")
