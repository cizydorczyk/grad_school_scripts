This script calculates the total bases in a fastq file by parsing the fastqc
data file (fastqc_data.txt).

It multiplies the number of reads of each length by the length, and returns the
total in a tsv file, with the first column being the filename and the second being
the total number of bases.

The script takes two arguments:

  fastqc_data_file: the fastqc_data.txt file for an isolate

  outputtsv: output .tsv file

  The output file does not have to exist when first running this script. That said,
  if the output file already does exist, the script will simply append to it if the
  same output file is used for several isolates.

**********Doesn't work if lengths are binned****************
