#function_1 extracts the filename and number of ottal seqeunce reads for each file it is given
#lists a and b are used to generate the list of fastqc_data.txt files, which
#will be iterated through after and function_1 will be called on each file
#the output is a single file, which is specified inside the function_1 defn

def function_1(filename1):
    f = open(filename1, 'r')
    for line in f:
        if line.startswith("Filename"):
            a = line.strip().split('\t')
        elif line.startswith("Total Sequences"):
            b = line.strip().split('\t')
    f.close()
    total = a + b
#Creates a new file in the specified directory and writes the filename and total reads to it as a .tsv file
    with open("/Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_sequence_reads/sequencereads.tsv", 'a') as new_file:
        new_file.write(total[1]+"\t")
        new_file.write(total[3]+"\n")


list1 = function_1("/Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_output_unzipped/6_R1_fastqc/fastqc_data.txt")

g = open("/Users/conradizydorczyk/grad_school_data/primary_project/file_list.txt", 'r')

#a is a list of all file numbers, used to create a list of all file paths for list b:
a = []
for line in g:
    a.append(line.strip())
g.close()

#b is a list of all unzipped fastqc data files, which will be iterated through
#and with each iteration, function_1 will extract the filename and # of reads
#and write to a singel file
b = []
for num in a:
    b.append("/Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_output_unzipped/"+num+"_R1_fastqc/fastqc_data.txt")
    b.append("/Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_output_unzipped/"+num+"_R2_fastqc/fastqc_data.txt")

#calling function_1 here on the file list stored in 'b'; output is a single tsv file
for i in b:
    function_1(i)
