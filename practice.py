from sys import argv
import os.path

script, krakensummary, output = argv
isolate = krakensummary.strip().split('_')[0]

class KrakenLine(object):
	def __init__(self, percent, numrootreads, numtaxonreads, rankcode, ncbitaxid, sciname, record):
		self.percent = percent
		self.numrootreads = numrootreads
		self.numtaxonreads = numtaxonreads
		self.rankcode = rankcode
		self.ncbitaxid = ncbitaxid
		self.sciname = sciname
		self.record = record

krakenline_objects = []
with open(krakensummary, 'r') as infile1:
	for line in infile1:
		temp1 = line.strip().split('\t')
		temp2 = temp1[5].strip()

		krakenline_objects.append(KrakenLine(float(temp1[0]), int(temp1[1]), int(temp1[2]), temp1[3], int(temp1[4]), temp2, temp1))

for i in krakenline_objects:
	if i.sciname == 'Pseudomonas aeruginosa':
		perc_pae = i.percent
	elif i.sciname == 'unclassified':
		perc_unclassified = i.percent

perc_other = 100 - perc_pae - perc_unclassified

print 'Percent Pae: ', perc_pae
print 'Percent Other: ', perc_other
print 'Percent Unclassified: ', perc_unclassified

if not os.path.isfile(output):
	with open(output, 'a+') as outfile1:
		outfile1.write('Isolate' + '\t' + 'Percent_Pae' + '\t' + 'Percent_Other' + '\t' + 'Percent_Unclass.' + '\n')
		outfile1.write(isolate + '\t' + str(perc_pae) + '\t' + str(perc_other) + '\t' + str(perc_unclassified) + '\n')
elif os.path.isfile(output):
	with open(output, 'a+') as outfile1:
		outfile1.write(isolate + '\t' + str(perc_pae) + '\t' + str(perc_other) + '\t' + str(perc_unclassified) + '\n')
