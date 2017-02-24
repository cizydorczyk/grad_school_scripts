"""This script takes two arguments: file with SNPs (.txt), and gene annotation file (tab separated)"""

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PSEUDOCODE:

# Create empty list of SNP positions
# Open file with SNPs:
#    For line in file:
#        If the line is not the header (first line) and not a duplicate:
#            Add SNP position to empty list of SNP positions

# Create empty list of gene annotation items
# Open file with gene annotations:
#    For line in file:
#        If line is not header (first line):
#            Create gene annotation item for line
#            Add item to empty list of annotation items

# For every SNP position:
#    For every gene annotation item:
#        If SNP position falls within gene:
#            Print gene description
#        Else:
#            Test next gene annotation item
#    If SNP position not found in any genes:
#        Print SNP is intergenic

# CODE BELOW:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# From sys module, import function argv
from sys import argv

# Define command line arguments (input files):
script, snp_file, annotation_file = argv

# Define function to obtain SNP positions from file with SNPs :
def identify_snp_positions(file_with_snps):
#   Create empty list for SNP positions:
    snp_positions_list = []
#   Open file with SNPs:
    with open(file_with_snps, 'r') as infile:
#       For every SNP in the SNP file, if it is not a duplicate, append its position to the empty list:
        for line in infile:
            if not line.startswith('Position') and line.split('\t')[0] not in snp_positions_list:
                snp_positions_list.append(int(line.split('\t')[0]))
#   Return list of SNP positions
    return snp_positions_list

snp_positions = identify_snp_positions(snp_file)

# Create class to define gene annotation objects:
class GeneAnnotation(object):
#   Initialize class with an argument for each piece of information in a gene annotation:
    def __init__(self, start_position, end_position, description, strand):
        self.start_position = start_position
        self.end_position = end_position
        self.description = description
        self.strand = strand

# Define function to create gene annotation class objects from gene annotation file:
def create_annotations(annotation_file):
#   Create empty list for gene annotation class objects:
    annotation_objects_list = []
#   Open gene annotation file:
    with open(annotation_file, 'r') as infile2:
#       For every gene annotation, create an instance of the GeneAnnotation class and append it to the empty list:
        for line in infile2:
            if not line.startswith('Start'):
                temp1 = line.strip('\n').split('\t')
                annotation_objects_list.append(GeneAnnotation(int(temp1[0]), int(temp1[1]), temp1[2], temp1[3]))
#   Return list of gene annotation objects
    return annotation_objects_list

annotations_list = create_annotations(annotation_file)

# Define function to test which gene a SNP occurs in, or whether it is intergenic:
def test_intergenic(list_of_snp_positions, list_of_GeneAnnotation_objects):
#   Test whether each SNP position falls within a gene:
    for position in list_of_snp_positions:
        for annotation in list_of_GeneAnnotation_objects:
#           If SNP position falls within gene, print gene description:
            if position in range(annotation.start_position, annotation.end_position+1):
                print "SNP at position " + str(position) + " is in gene: " + str(annotation.description)
                break
#           If SNP does not fall within gene, test next gene:
            else:
                continue
#       If SNP does not fall within any genes, print SNP is intergenic:
        else:
            print "SNP at position " + str(position) + " is intergenic."

test_intergenic(snp_positions, annotations_list)
