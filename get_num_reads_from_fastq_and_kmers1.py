from sys import argv

# outputtsv must already exist! Use bash (touch) to make!
script, fastqfile, kmerfile, outputtsv = argv

fastqcount = 0
with open(fastqfile, 'r') as fastqfile1:
    for line in fastqfile1:
        if line.startswith("@"):
            fastqcount += 1
        else:
            pass

kmercount = 0
with open(kmerfile, 'r') as kmerfile:
    for line in kmerfile:
        if line.startswith("@"):
            kmercount += 1
        else:
            pass

with open(outputtsv, 'a') as outputfile:
    outputfile.write(fastqfile + "\t")
    outputfile.write(str(fastqcount)+"\t")
    outputfile.write(str(kmercount)+"\n")
