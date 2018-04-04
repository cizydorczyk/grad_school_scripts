from sys import argv
import re
from ete3 import Tree

script, inputtree = argv

t = Tree('((((H,K)D,(F,I)G)B,E)A,((L,(N,Q)O)J,(P,S)M)C);', format=1)

for node in t.traverse("preorder"):
    print node


with open(inputtree, 'r') as infile1:
    for line in infile1:
        line_list = line.strip().split(',')

line_list2 = []
for i in line_list:
    i = re.sub("\(", "", i)
    i = re.sub("\)", "", i)
    i = re.sub(";", "", i)
    i_list = i.split(':')
    for j in i_list:
        if "pae" not in j:
            line_list2.append(float(j))

total_branch_length = 0
for length in line_list2:
    total_branch_length += length

print total_branch_length
print 1/total_branch_length