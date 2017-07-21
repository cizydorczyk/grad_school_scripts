from sys import argv

script, inputblast, inputfasta, outputfasta = argv

contig_list = []
with open(inputblast, 'r') as infile:
    for line in infile:
        if "Filename" not in line:
            line_ = line.strip().split('\t')
            if 'orange' not in line_:
                if float(line_[2]) > 10.0:
                    contig_list.append(line_[0])

contig_seq_dict = {}
with open(inputfasta, 'r') as infile2:
    fasta = infile2.read()
    fasta = fasta.split(">")
    for i in fasta:
        new_line = i.split("\n")
        header = new_line[0]
        new_line.pop(0)
        sequence = ''.join(new_line)
        contig_seq_dict[header] = sequence

with open(outputfasta, 'w') as outfile:
    for key, value in contig_seq_dict.iteritems():
        if key in contig_list:
            outfile.write(">" + key + '\n' + value + '\n')
        else:
            pass
