from sys import argv

# inputblasefile is the input blast output file (tsv)
# inputcontigfasta is the fasta file containing all contigs for an isolate
# outputfafile is the output fasta file containing only contigs that did not
#   blast to Pseudomonas aeruginosa
# outputcontfile is a tsv file containing the contig id and organism of contigs
#   that did not blast to Pseudomonas aeruginosa
script, inputblastfile, inputcontigfasta, outputfafile, outputcontfile = argv

#print inputblastfile
#print inputcontigfasta

print "Now processing "+inputblastfile

inputblast = open(inputblastfile, 'r')

contfile = open(outputcontfile, 'w')

contigs_to_remove = []
for line in inputblast:
    if not line.startswith("#"):
        temp1 = line.strip("\n").split("\t")
        if not temp1[-1].startswith("Pseudomonas aeruginosa"):
            contigs_to_remove.append(">"+temp1[0])
            contfile.write(temp1[0]+"\t")
            contfile.write(temp1[-1]+"\n")
inputblast.close()
contfile.close()

outfasta = open(outputfafile, 'w')

contigs = [line.strip("\n") for line in open(inputcontigfasta, 'r')]

infile = open(inputcontigfasta, 'r')

for line in infile:
    if line.strip("\n") in contigs_to_remove:
        outfasta.write(line)
        line = next(infile)
        while not line.startswith(">") and line != '':
            outfasta.write(line)
            line = next(infile)

infile.close()


#outfasta.close()



#contigfile = open(inputcontigfasta, 'r')
#contigs = []
#for line in contigfile:
#    temp2 = line.strip("\n")
#    contigs.append(temp2)
#contigfile.close()

#for inum, i in enumerate(contigs):
#    if i in contigs_to_remove:
#        count = 1
#        outfile.write(i+"\n")
#        while not contigs[inum+count].startswith(">"):
#            outfile.write(contigs[inum+count]+"\n")
#            count += 1

#outfile.close()
