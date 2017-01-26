This script takes in two corresponding fastq files (R1 and R2) for an isolate and removes
reads from them that are found in two corresponding contaminant/kmer fastq files.

It first parses both R1 and R2 fastq files, converting each read into a read object, and
creating a list of read objects for each file. Each list is then converted to an ordered set
that maintains the order of read objects (items).

It then iterates through both contaminant/kmer fastq files, and creates a single, combined
list of headers from both. These are the reads that will be removed from the two input fastq
files. It then converts this single list into an ordered set as well (though order isn't
important here except for speed perhaps...).

It then iterates through the two intput fastq file sets (for R1 and R2) simultaneously.
If the header of R1 for each iteration is not in the list of headers to remove, it
appends the full read (R1) to a list, and the full read (R2) to a separate list.
After iterating through all read objects in R1, it writes the resulting lists to
two separate output files.

The script takes six argumnets as inputs:

    fastq1 = R1 fastq file for isolate (same isolate as fastq2)
    fastq2 = R2 fastq file for isolate (same isolate as fastq1)
    kmers1 = fastq file with reads with kmers for R1
    kmers2 = fastq file with reads with kmers for R2
    output1 = R1 output fastq file
    output2 = R2 output fastq file

None of the output files have to exist before running the script, and are overwritten
if they do exist.
