from sys import argv
from collections import OrderedDict
import pandas as pd
from scipy.stats import fisher_exact

script, input_alignment_file, input_positions_file, output_file, output_contingency_tables = argv

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

# Create list of positions in reference, and reference sequence:
positions_list = []
reference = ''
with open(input_positions_file, 'r') as infile2:
    for line in infile2:
        if not line.startswith("Position"):
			reference += line.strip().split('\t')[1]
			positions_list.append(int(line.strip().split('\t')[0]))

# convert snp sequences to 0 (reference base) and 1 (alt base):
snp_matrix_dict = {}
for i in fasta_seq:
	sequence = ''
	for j in range(0,len(reference)):
		if reference[j] == fasta_seq[i][j]:
			sequence += '0'
		else:
			sequence += '1'
	snp_matrix_dict[i] = list(sequence)

# Turn sequence dictionary into dataframe, with high quality positions as the column names:
df1 = pd.DataFrame.from_dict(snp_matrix_dict, orient='index')
df1.columns = positions_list

df1['Phenotype'] = pd.Series(eradication_dict)

# with open(output_file, 'w') as outfile1:
output_list = []
output_contingency = []

for i in positions_list:
	tab = pd.crosstab(df1[i] == '1', df1["Phenotype"] == '1')
	output_contingency.append(tab.to_string())
	print tab
	try:
		oddsratio, p_value = fisher_exact(tab)
		output_list.append(str(i) + '\t' + str(oddsratio) + '\t' + str(p_value))
	except ValueError:
		oddsratio = "NA"
		p_value = "NA"
		output_list.append(str(i) + '\t' + oddsratio + '\t' + p_value)

with open(output_file, 'w') as outfile1:
	outfile1.write('SNP' + '\t' + 'oddsratio' + '\t' + 'p_value' + '\n')
	outfile1.write('\n'.join(output_list))

with open(output_contingency_tables, 'w') as outfile2:
	outfile2.write('\n'.join(output_contingency))
