This script takes a quast report.tsv file and parses it, counting how many
assembly parameters are better in assembly 1 vs assembly 2. This script
was designed initially to look at qc5 vs qc20 trimmed assemblies. It should
be run as part of a loop for ALL quast report.tsv files, as the output is a
single tsv file with the number of parameters that were better for each
assembly.

This script takes 3 arguments:

  isolate_number: "isolate_"$i, where $i should be the isolate # from the loop
  inputfile1: input report.txt file for an isolate
  outputtsv: the output tsv file

If the output tsv file does not exist prior to running the script, it is
created. If it exists, it is simply appended to.
