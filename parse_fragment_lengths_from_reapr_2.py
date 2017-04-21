from sys import argv
import os.path
import numpy
import math


script, input_fragment_lengths, output = argv

fragment_lengths = []
fragment_counts = []

with open(input_fragment_lengths, 'r') as infile:
    for line in infile:
        fragment_lengths.append(int(line.strip().split('\t')[0]))
        fragment_counts.append(int(line.strip().split('\t')[1]))

lengths_array = numpy.array(fragment_lengths)
counts_array = numpy.array(fragment_counts)

weighted_average = numpy.average(a=lengths_array, weights=counts_array)
weighted_stdev = math.sqrt(numpy.average((lengths_array-weighted_average)**2, weights=counts_array))

if os.path.isfile(output):
    with open(output, 'a') as outfile:

        outfile.write(input_fragment_lengths.strip().split("/")[-2] + '\t' + str(weighted_average) + '\t' + str(weighted_stdev) + '\n')

elif not os.path.isfile(output):
    with open(output, 'w') as outfile:
        outfile.write("Isolate" + '\t' + "Mean" + '\t' + "Stdev" + '\n')
        outfile.write(input_fragment_lengths.strip().split("/")[-2] + '\t' + str(weighted_average) + '\t' + str(weighted_stdev) + '\n')
