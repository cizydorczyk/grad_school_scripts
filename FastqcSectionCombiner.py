#!/usr/bin/env python

import sys, getopt
import os.path

HELP="""
Program: FastqcSectionCombiner
Version: 0.1
Contact: Julio Diaz Caballero <julio.diaz@mail.utoronto.ca>
Details: Combines sections obtained in FastqcSectionGetter from various
         input files. The sections should be the same.

Usage:\tFastqcSectionCombiner [-v|-h] -i <input> -o <output.tsv>

Arguments:\t-i\tInput file including names of section files from
\t\t\tdifferent input files created with FastqSectionGetter
\t\t-o\tOutput file

Options:\t-v\tverbose output
\t\t-h\thelp message\n"""

inputfile = ''
outputfile = ''
argv=sys.argv[1:]
verbose=False
if(len(argv)==0):
    print HELP
    sys.exit(1)
try:
    opts, args = getopt.getopt(argv,"hvi:o:",["ifile=","ofile="])
except getopt.GetoptError:
    print HELP
    sys.exit(1)
for opt, arg in opts:
    if opt == '-h':
        print HELP
        sys.exit(0)
    elif opt == '-v':
        verbose=True
    elif opt in ("-i", "--ifile"):
        inputfile = arg
    elif opt in ("-o", "--ofile"):
        outputfile = arg

if(inputfile==''):
    print "\n[Needs an input file]\n"
    sys.exit(1)
elif(os.path.isfile(inputfile)==False):
    print "\n[Input file does not exist]\n"
    sys.exit(1)
if(outputfile==''):
    print "\n[Needs an output file]\n"
    sys.exit(1)

if(os.path.isdir(os.path.dirname(outputfile))==False and os.path.dirname(outputfile)!=""):
    print "\n[Directory does not exist]\n"
    sys.exit(1)
if(os.path.isdir(outputfile)):
    print "\n[File is already a directory]\n"
    sys.exit(1)


fileList = list()
with open(inputfile ,"r") as lns:
    for fline in lns:
        fileList.append(fline.rstrip())

table = {}

nlist = list()
with open(fileList[0], "r") as nme:
    for nline in nme:
        if(nline[0]!='#'):
            nlist.append(nline.rstrip().split()[0])
table[""] = nlist

for fileName in fileList:
    my_list = list()
    with open(fileName, "r") as ins:
        for line in ins:
            if(line[0]!='#'):
                my_list.append(line.rstrip().split()[1])
    table[fileName] = my_list

#Print as table
f= open(outputfile, 'w')

for key in table.keys():
    f.write(key+"\t")
f.write("\n")
#    sys.stdout.write(key+"\t")
#sys.stdout.write("\n")

for num in range(0,len(table[table.keys()[0]])):
    for key in table.keys():
        f.write(table[key][num]+"\t")
    f.write("\n")
#        sys.stdout.write(table[key][num]+"\t")
#    sys.stdout.write("\n")

f.close()
