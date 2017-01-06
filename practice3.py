

import sys

class Variants(object):
    def __init__(self, CHROM, POS, REF, ALT, QUAL, DP, DP4):
        self.CHROM = CHROM
        self.POS = POS
        self.REF = REF
        self.ALT = ALT
        self.QUAL = QUAL
        #self.INFO = INFO
        self.DP = DP
        self.DP4 = DP4

    def print_all(self):
        print self.CHROM + ' ' + str(self.POS) + ' ' + self.REF + ' ' + self.ALT + ' ' + str(self.QUAL) + ' ' + str(self.DP) + ' ' + str(self.DP4)


#var_lst is a list of ALL variants w/ specific info I want
var_lst = []

#variant_list is a list of Variant class objects, one for each item in var_lst
variant_list = []

for i in sys.argv:
    if i != "practice3.py":

        f = open(i, "r")

        for line in f:
            if not "#" in line:
                print line
        #        #new_lst contains only elements I want from a variant line
        #        new_lst = []
        #        new_lst.append(temp_lst[0])
        #        new_lst.append(int(temp_lst[1]))
        #        new_lst.append(temp_lst[3])
        #        new_lst.append(temp_lst[4])
        #        new_lst.append(float(temp_lst[5]))
        #        temp2_lst = temp_lst[7].split(";")
        #        #append just the numerical DP value:
        #        for i in temp2_lst:
        #            if "DP=" in i:
        #                if len(i) >= 5:
        #                    new_lst.append(int(temp2_lst[temp2_lst.index(i)][-2] + temp2_lst[temp2_lst.index(i)][-1]))
        #                elif len(i) == 4:
        #                    new_lst.append(int(temp2_lst[temp2_lst.index(i)][-1]))
                            #append int DP4 values in the form of a list to new_lst:
        #        for i in temp2_lst:
        #             if "DP4" in i:
        #                 temp3 = i[4:].split(",")
        #                 temp4 = []
        #                 for i in temp3:
        #                     temp4.append(int(i))
        #                 new_lst.append(temp4)

        #        var_lst.append(new_lst)


        f.close()

#create instance of Variant class for each line (item in var_lst)
#and add it to variant_list:
#for i in var_lst:
#    temp1 = Variants(i[0], i[1], i[2], i[3], i[4], i[5], i[6])
#    variant_list.append(temp1)

#indel_list contains a list of positions with indels:
#indel_list = []
#add any indel positions to indel_list:
#for i in variant_list:
#    if len(i.REF) > 1 or len(i.ALT) > 1:
#        indel_list.append(int(i.POS))

#filtered_variant_POS_list is a list of positions filtered with the following
#parameters: DP > 20; QUAL > 30; REF forward+reverse is less than 15% of
#ALT forward+reverse; ALT forward and reverse are at least 3; at least 150bp
#away from ends of chromosome (reference genome; 6264044bp Genbank AE004091.2)
#or any indels
#filtered_variant_POS_list = []

#for i in variant_list:
#    if i.QUAL >= 30 and i.DP >= 20 and i.DP4[2] >= 3 and i.DP4[3] >= 3 and (i.DP4[0] + i.DP4[1]) < 0.15*(i.DP4[2] + i.DP4[3]) and i.POS > 150 and i.POS < (6264044 - 150):
#        for j in indel_list:
#            if i.POS not in range(j-150, j+151) and i.POS not in filtered_variant_POS_list:
#                filtered_variant_POS_list.append(i.POS)

#print filtered_variant_POS_list
