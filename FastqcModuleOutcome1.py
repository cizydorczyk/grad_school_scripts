g = open("/Users/conradizydorczyk/grad_school_data/primary_project/file_list.txt", 'r')

#a is a list of all file numbers:
a = []
for line in g:
    a.append(line.strip())
g.close()

#b is a list of all unzipped fastqc data files:
b = []
for num in a:
    b.append("/Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_output_unzipped/"+num+"_R1_fastqc/fastqc_data.txt")
    b.append("/Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_output_unzipped/"+num+"_R2_fastqc/fastqc_data.txt")

def fastqc_parser1(filelist):
    header = ["Filename", "Per base sequence quality", "Per sequence quality scores", "Per base sequence content",
              "Per sequence GC content", "Per base N content", "Sequence Length Distribution",
              "Sequence Duplication Levels", "Overrepresented sequences", "Adapter Content", "Kmer Content"]
#make sure the file open here is in the directory you want and make sure you open this same file again later on line 62
    with open("/Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_module_outcomes.tsv", 'w') as testfile:
        for i in header[0:-1]:
            testfile.write(i+'\t')
        testfile.write(header[-1]+'\n')
    for item in filelist:
        module_outcome = []
        f = open(item, 'r')
        for line in f:
            if line.startswith("Filename"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Per base sequence quality"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Per sequence quality scores"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Per base sequence content"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Per sequence GC content"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Per base N content"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Sequence Length Distribution"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Sequence Duplication Levels"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Overrepresented sequences"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Adapter Content"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
            elif line.startswith(">>Kmer Content"):
                a = line.strip().split('\t')
                module_outcome.append(a[-1])
#make sure the file open here matches the one open on line 20
        with open("/Users/conradizydorczyk/grad_school_data/primary_project/Fastqc/Fastqc_raw/fastqc_module_outcomes.tsv", 'a') as testfile2:
            for j in module_outcome[0:-1]:
                testfile2.write(j+'\t')
            testfile2.write(module_outcome[-1]+'\n')

fastqc_parser1(b)