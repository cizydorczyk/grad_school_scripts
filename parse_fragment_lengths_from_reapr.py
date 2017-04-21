from sys import argv
import os.path

script, input_fragment_lengths, output = argv

dist = {}

with open(input_fragment_lengths, 'r') as infile:
    for line in infile:
        dist[(int(line.strip().split('\t')[0]))] = (int(line.strip().split('\t')[1]))

bins = {"0-49":0, "050-99":0, '100-149':0, '150-199':0, '200-249':0, '250-299':0, '300-349':0, '350-399':0, '400-449':0, '450-499':0, '500-549':0, '550-599':0, '600-649':0, '650-699':0, '700-749':0, '750-799':0, '750-799':0, '800-849':0, '850-899':0, '900-949':0, '950-999':0, '>= 1000':0}


for key, value in dist.iteritems():
    if key <= 49:
        bins["0-49"] += value
    elif 50 <= key <= 99:
        bins['050-99'] += value
    elif 100 <= key <= 149:
        bins['100-149'] += value
    elif 150 <= key <= 199:
        bins['150-199'] += value
    elif 200 <= key <= 249:
        bins['200-249'] += value
    elif 250 <= key <= 299:
        bins['250-299'] += value
    elif 300 <= key <= 349:
        bins['300-349'] += value
    elif 350 <= key <= 399:
        bins['350-399'] += value
    elif 400 <= key <= 449:
        bins['400-449'] += value
    elif 450 <= key <= 499:
        bins['450-499'] += value
    elif 500 <= key <= 549:
        bins['500-549'] += value
    elif 550 <= key <= 599:
        bins['550-599'] += value
    elif 600 <= key <= 649:
        bins['600-649'] += value
    elif 650 <= key <= 699:
        bins['650-699'] += value
    elif 700 <= key <= 749:
        bins['700-749'] += value
    elif 750 <= key <= 799:
        bins['750-799'] += value
    elif 800 <= key <= 849:
        bins['800-849'] += value
    elif 850 <= key <= 899:
        bins['850-899'] += value
    elif 900 <= key <= 949:
        bins['900-949'] += value
    elif 950 <= key <= 999:
        bins['950-999'] += value
    elif 1000 <= key:
        bins['>= 1000'] += value

header = []
isolate_no = "isolate_number"
isolate_num = input_fragment_lengths.split('/')[-2]

header.append(isolate_no)

for i in sorted(bins):

    header.append(i)

if os.path.isfile(output):
    with open(output, 'a') as outfile:
        values = []
        values.append(isolate_num)
        for i in sorted(bins):
            values.append(str(bins[i]))
        outfile.write('\t'.join(values) + '\n')
elif not os.path.isfile(output):
    with open(output, 'w') as outfile:
        outfile.write('\t'.join(header) + '\n')
        values = []
        values.append(isolate_num)
        for i in sorted(bins):
            values.append(str(bins[i]))
        outfile.write('\t'.join(values) + '\n')
