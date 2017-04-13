from sys import argv

script, inputfile, newoutputfile = argv
#print inputfile,

list1 = []

with open(inputfile, 'r') as infile:
    for line in infile:
        if "Filename" not in line:
            list1.append(line.strip())

with open(newoutputfile, 'w') as outfile:
    outfile.write("Filename" + '\t' + "Length" + '\t' + "Coverage" + '\t' + "Color" + '\t' + "Outline" + '\n')
    for i in list1:
        outfile.write(i + '\n')
