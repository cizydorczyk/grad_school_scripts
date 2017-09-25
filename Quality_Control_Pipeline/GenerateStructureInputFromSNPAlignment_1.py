from sys import argv
from collections import OrderedDict
import pandas as pd
from scipy.stats import fisher_exact

script, input_alignment_file, input_positions_file, output_file = argv

# Input positions file is a tab separated file with position\tref_base format, no header

eradication_dict = {'6':'0', '7':'0', '8':'0', '10':'0', '14':'0', '17':'0', '19':'0',
'21':'0', '23':'0', '30':'0', '36':'0', '37':'0', '53':'0', '97':'0', '130':'0', '152':'0', '160':'0',
'180':'0', '196':'0', '200':'0', '217':'0', '225':'0', '235':'0', '251':'0', '257':'1',
'258':'1', '259':'0', '263':'0', '265':'0', '270':'0', '273':'0', '276':'0', '282':'0',
'275':'0', '278':'0', '280':'0', '288':'0', '290':'0', '293':'0', '298C':'1', '303':'0',
'306':'0', '310':'0', '315':'0', '317':'0', '318':'0', '323':'0', '325':'0', '326':'0',
'327':'0', '330':'1', '331':'1', '332':'1', '335':'1', '336':'1', '337':'1', '340':'1',
'341':'1', '342':'1', '358':'0', '359':'0', '360':'0', '366':'0', '367':'0', '369':'0',
'375':'1', '379':'1', '380':'1', '385':'1', '388':'1', '390':'0', '392':'0', '393':'0',
'395':'0', '400':'0', '401':'0', '402':'0', '404C':'0', '405C':'0', '406C':'0', '409':'0',
'410C':'0', '411C':'0', '412':'0', '417C':'0', '418C':'0', '419C':'0', '420N2':'0',
'427':'0', '428':'0', '429':'0', '430':'0', '431':'0', '432':'0', '433':'0', '446':'0',
'448':'0', '449':'0', '450':'0', '453':'0', '455':'0', '457':'0', '459C':'0', '465':'1',
'471':'0', '472':'0', '473':'0', '475':'0', '476':'0', '479':'0', '480':'0', '487':'0',
'501':'0', '503':'1', '504':'1', '505':'1', '506':'1', '507':'0', '508':'0', '509':'0',
'510':'0', '511':'0', '512':'0', '513':'0', '514':'0', '521':'0', '525':'0', '526':'0',
'527':'0', '539':'0', '540':'0', '541-2':'0', '544':'0', '551':'1', '552':'1', '553C':'1',
'557':'0', '559-2':'0', '562':'0', '563':'0', '566':'1', '567':'1', '568':'1', '569':'0',
'570':'0', '571':'0', '572-1':'1', '573-2':'1', '573-3':'1', '575':'0', '577':'0', '578':'0',
'580':'1', '581':'1', '583':'1'}

print len(eradication_dict)
count = 0
for i in eradication_dict:
	if eradication_dict[i] == '1':
		count += 1
print count

fasta_seq ={}

with open(input_alignment_file, 'r') as infile1:
    for line in infile1:
        if line.startswith(">"):
            fasta_seq[line.strip()[1:]] = list(next(infile1).strip())

# Create list of high quality positions:
positions_list = []
positions_str_list = []

with open(input_positions_file, 'r') as infile2:
    for line in infile2:
        if not line.startswith("Position"):
            positions_list.append(int(line.strip().split('\t')[0]))
            positions_str_list.append(line.strip().split('\t'[0]))

# Turn sequence dictionary into dataframe, with high quality positions as the column names:
df1 = pd.DataFrame.from_dict(fasta_seq, orient='index')
df1.columns = positions_list

print df1[0:10]

df1.replace(['A', 'C', 'G', 'T'], [100, 200, 300, 400], inplace=True)


#df1.to_csv(path_or_buf=output_file, sep='\t')
df2 = df1.ix[:,[154,155]]
df2.to_csv(path_or_buf=output_file, sep='\t')
# Convert dataframe into dictionary, with column name as key and list of values in each column as value:
# columns = df1.to_dict(orient='list')



# columns2 = {}
# for i in sorted(columns):
#     iset = set(columns[i])
#     isetlist = list(iset)
#     new_list = []
#     if len(isetlist) == 2:
#         # print "lenset 2"
#         for j in columns[i]:
#             if isetlist[0] == j:
#                 new_list.append(0)
#             else:
#                 new_list.append(1)
#     elif len(isetlist) == 3:
#         # print "lenset 3"
#         for j in columns[i]:
#             if isetlist[0] == j:
#                 new_list.append(0)
#             elif isetlist[1] == 1:
#                 new_list.append(1)
#             else:
#                 new_list.append(2)
#     elif len(isetlist) == 4:
#         # print "lenset 4"
#         for j in columns[i]:
#             if isetlist[0] == j:
#                 new_list.append(0)
#             elif isetlist[1] == j:
#                 new_list.append(1)
#             elif isetlist[2] == j:
#                 new_list.append(2)
#             else:
#                 new_list.append(3)
#
#     columns2[i] = new_list
#
# for i in sorted(columns2):
#     print i, columns2[i]

# with open(output_file, 'w') as outfile1:
#     outfile1.write('\t' + '\t'.join(positions_str_list) + '\n')
#     to_write = []
#     for i in sorted(columns2):
#         to_write.append(str(i) + '\t' + '\t'.join(columns2[i]) + '\n')
#     outfile1.write(to_write)
