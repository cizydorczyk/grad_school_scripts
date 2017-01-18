from sys import argv
import itertools
import pandas

script, inputblastfile, inputcontigstats, output_tsv = argv


def blast_parse(inputblastfile):
    print "Working on file: " + str(inputblastfile)
    data = pandas.read_table(inputcontigstats, sep="\t")
    length_list = list(data["Consensus length"])
    coverage_list = list(data["Average coverage"])
    with open(inputblastfile, 'r') as infile:
        lines = list(infile)
    query_list = []
    top_hit_list = []
    fill_color = []
    out_color = []
    for lnum, i in enumerate(lines):
        if "# Query: " in i:
            if (lnum+4) < len(lines):
                temp1 = i.split(' ')
                query_list.append(temp1[-1].strip("\n"))
                print "\tQuery at index: " + str(lnum)
                temp2 = lines[lnum+4]
                temp3 = temp2.split("\t")
                top_hit_list.append(temp3[-1].strip("\n"))
                if "Pseudomonas aeruginosa" in temp3[-1]:
                    fill_color.append("deepskyblue")
                    out_color.append("blue")
                else:
                    fill_color.append("yellow")
                    out_color.append("orange")
            elif (lnum+4) >= len(lines):
                print "\tQuery at final index: " + str(lnum)
                temp1 = i.split(' ')
                query_list.append(temp1[-1].strip("\n"))
                fill_color.append("yellow")
                out_color.append("orange")
                top_hit_list.append("N/A")
    with open(output_tsv, 'w') as outfile:
        for query, length, coverage, fill, out, hit in itertools.izip(query_list, length_list, coverage_list, fill_color, out_color, top_hit_list):
            outfile.write(str(query) + "\t" + str(length) + "\t" + str(coverage) + "\t" + str(fill) + "\t" + str(out) + "\t" + str(hit) + "\n")




contigs_to_remove = blast_parse(inputblastfile)
