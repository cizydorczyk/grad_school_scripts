from sys import argv

script, inputfile = argv

file_list = []
with open(inputfile, 'r') as infile1:
    for line in infile1:
        file_list.append(line.strip())

req_elements = ['d', 'p', 'c', 'o', 'f', 'g']

output_list = []

for record in file_list:
    record_list = record.split(',')
    
    if len(record_list) == 7:
        # print ','.join(record_list)
        output_list.append(','.join(record_list))

    elif len(record_list) < 7:
        record_elements = []

        for i in record_list[1:]:
            record_elements.append(i.split('__')[0])
        print record_elements
        
        # missing_elements = [i for i in req_elements if i not in record_elements]
        
        # for elem in missing_elements:
        #     if elem == 'd':
        #         record_list.insert(1, 'NA')
        #     elif elem == 'p':
        #         record_list.insert(2, 'NA')
        #     elif elem == 'c':
        #         record_list.insert(3, 'NA')
        #     elif elem == 'o':
        #         record_list.insert(4, 'NA')
        #     elif elem == 'f':
        #         record_list.insert(5, 'NA')
        #     elif elem == 'g':
        #         record_list.insert(6, 'NA')
        
        # # print ','.join(record_list)
        # output_list.append(','.join(record_list))

print '\n'.join(output_list)