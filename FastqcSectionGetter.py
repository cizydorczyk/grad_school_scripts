#!/usr/bin/env python

#from sys import argv
import sys, getopt
import os.path

HELP="""
Program: FastqcSectionGetter
Version: 0.1
Contact: Julio Diaz Caballero <julio.diaz@mail.utoronto.ca>
Details: Gets a specific section of the fastqc report.

Usage:\tFastqcSectionGetter [-v|-h] <read.fastq> [-o <output.fastq>]

Arguments:\t-i\tInput file. The input file should be 'fastqc_data.txt'
\t\t\tas created by fastqc
\t\t-s\tSection of fatqc to extract:
\t\t\tadapter\t\tAdapter content
\t\t\tquality\t\tPer base sequence quality
\t\t\tbase_content\tPer base sequence content
\t\t\tgc_content\tPer sequence GC content
\t\t\tn_count\t\tPer base N content
\t\t\tlength_dist\tSequence length distribution
\t\t\tduplicate\tSequence duplication levels
\t\t\tkmer\t\tKmer content

Options:\t-o\tOutput file. If this option is not set, the output
\t\t\tfile will be named according to the section
\t\t-v\tverbose output
\t\t-h\thelp message\n"""

LENGTH=151

inputfile = ''
outputfile = ''
section = ''
sectionList={"adapter":"Adapter Content",
             "quality":"Per base sequence quality",
             "base_content":"Per base sequence content",
             "gc_content":"Per sequence GC content",
             "n_count":"Per base N content",
             "length_dist":"Sequence Length Distribution",
             "duplicate":"Sequence Duplication Levels",
             "kmer":"Kmer Content"}
argv=sys.argv[1:]
verbose=False
if(len(argv)==0):
    print HELP
    sys.exit(1)
try:
    opts, args = getopt.getopt(argv,"hvs:i:o:",["ifile=","ofile="])
except getopt.GetoptError:
    print HELP
    sys.exit(1)
for opt, arg in opts:
    if opt == '-h':
        print HELP
        sys.exit(0)
    elif opt == '-v':
        verbose=True
    elif opt == '-s':
        section = arg
    elif opt in ("-i", "--ifile"):
        inputfile = arg
    elif opt in ("-o", "--ofile"):
        outputfile = arg

if(section==''):
    print "\n[Needs to specify section]\n"
    sys.exit(1)
elif(section not in list(sectionList.keys())):
    print "\n[Unaccepted section]\tAccepted sections: \n"
    sys.exit(1)
if(inputfile==''):
    print "\n[Needs an input file]\n"
    sys.exit(1)
elif(os.path.isfile(inputfile)==False):
    print "\n[Input file does not exist]\n"
    sys.exit(1)
if(outputfile==''):
    outputfile=os.path.splitext(inputfile)[0]+"."+section+os.path.splitext(inputfile)[1]

if(os.path.isdir(os.path.dirname(outputfile))==False and os.path.dirname(outputfile)!=""):
    print "\n[Directory does not exist]\n"
    sys.exit(1)
if(os.path.isdir(outputfile)):
    print "\n[File is already a directory]\n"
    sys.exit(1)

if(verbose):
    print "\n[Input file found]"

f= open(outputfile, 'w')
module=sectionList[section]

if(verbose):
    print "[Searching for section]"

stillOn=False
with open(inputfile, "r") as ins:
    for line in ins:
        lineT=line.rstrip()
        if(lineT.split("\t")[0]==">>"+module):
            stillOn=True
            continue;
        if(stillOn==True and lineT!=">>END_MODULE"):
            #print lineT
            f.write(lineT+"\n")
        else:
            stillOn=False
f.close()

if(verbose):
    print "[Output written]\t"+section+" section was written in '"+os.path.basename(outputfile)+"'\n"
