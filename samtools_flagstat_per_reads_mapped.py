from sys import argv

"""Flagstat output file must have had "echo "new file\t"$i >> output.txt" line before each iteration of flagstat"""
script, flagstatoutput, outputfiletsv = argv

def parse(flagstatoutputfile, outputfiletsv):
    with open(flagstatoutputfile, 'r') as infile:
        lines = [line.strip("\n") for line in infile]
    count = 0
    with open(outputfiletsv, 'w') as outfile:
        for lnum, item in enumerate(lines):
            if item.startswith("new"):
                count += 1
                temp1 = item.split(" ")
                temp2 = lines[lnum+5].split(" ")
                temp3 = temp2[-3][1:-1]
                outfile.write(str(temp1[1])+"\t"+str(temp3)+"\n")
#                print temp1[1], temp3

    print "Files parsed: " + str(count)

parse(flagstatoutput, outputfiletsv)
