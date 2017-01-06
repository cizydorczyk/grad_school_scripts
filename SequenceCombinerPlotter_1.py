import csv

f = open("/home/conrad/Data/160902_NextSeq/22_eradicated_fastqc_output/FastqcSequenceGetterOutput/22_eradicated_total_quality.tsv", "rb")

csv_f = csv.reader(f, delimiter='\t')
three02 = [i[1] for i in csv_f]
print three02

f.close()
