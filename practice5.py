from sys import argv

script, output_file = argv

samples = range(1, 2251)

out_list = []
for number in samples:
    muscle_command = "/home/conrad/Software/muscle3.8.31_i86linux64 -in /home/conrad/Data/primary_project_3/ortholog_analysis/G66_ortho_families/core_fam_ffn_records/core_fam_" + str(number) + ".ffn -out /home/conrad/Data/primary_project_3/ortholog_analysis/G66_ortho_families/core_fam_ffn_alignments/core_fam_" + str(number) + ".fa"
    out_list.append(muscle_command)

with open(output_file, 'w') as outfile:
    outfile.write('\n'.join(out_list))
