from sys import argv
import itertools

script, inputblastfile, output_tsv = argv

#print inputblastfile
#print inputcontigfasta


def blast_parse(inputblastfile):
    with open(inputblastfile, 'r') as infile:
        lines = list(infile)
    query_list = []
    top_hit_list = []
    fill_color = []
    out_color = []
    for lnum, i in enumerate(lines):
        if "# Query: " in i:
            temp1 = i.split(' ')
            query_list.append(temp1[-1].strip("\n"))
            temp2 = lines[lnum+4]
            temp3 = temp2.split("\t")
            top_hit_list.append(temp3[-1].strip("\n"))
            if "Pseudomonas aeruginosa" in temp3[-1]:
                fill_color.append("deepskyblue")
                out_color.append("blue")
            else:
                fill_color.append("yellow")
                out_color.append("orange")
    with open(output_tsv, 'w') as outfile:
        for i, j, k, l in itertools.izip(query_list, top_hit_list, fill_color, out_color):
            outfile.write(str(i + "\t" + j + "\t" + k + "\t" + l + "\n"))

contigs_to_remove = blast_parse(inputblastfile)
