This script calls SNPS from variants identified by samtools/bcftools. It takes the following output:

	text file list of vcf files (full paths) filtered for high quality positions (filtered for minimum quality 30) for the bwa aligner
	text file list of vcf files (full paths) filtered for high quality positions (filtered for minimum quality 30) for the novoalign aligner
	text file list of vcf files (full paths) filtered for high quality positions (filtered for minimum quality 30) for the lastalign aligner

	Reference genome fasta file - must be the same as the reference fasta genome used for reference alignment
	An output file to write all high quality positions identified from above three aligners

	text file list of vcf files (full paths) for variants filtered with lower quality cutoff (minimum quality 25) for the bwa aligner
	text file list of vcf files (full paths) for variants filtered with lower quality cutoff (minimum quality 25) for the novoalign aligner
	text file list of vcf files (full paths) for variants filtered with lower quality cutoff (minimum quality 25) for the lastlaign aligner

	Output directory to which to write files containing list of SNPs identified for each isolate (all will be the same, but serves as a record of which positions
		were called as SNP, REF, and N for that isolate)

	Output fasta file that will contain the SNP alignment sequences, in fasta format

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
How the script works (pseudocode):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for each file list of quality 30 vcf files:
	open file and append each line to a separate list

Create master set of common positions

for each isolate:
  for each set of 3 quality 30 vcf files:
    identify common positions
    for each common position:
      if common position not in master set of common positions:
        append common position to master set

Create master dictionary for positions and reference bases at that position

With the reference fasta file open:
  Identify base at each position in reference sequence
  Append to master dictionary, with key = position and value = reference base

write master dictionary to file # to keep as a record of positions and ref bases

# Have a complete set of common positions now, where each position appears as a
# variant in at least one isolate, and all three aligners agree on that position
# being a variant, for that isolate

# Now to actually call SNPs for each isolate and generate alignment sequences:

For each isolate:
  Identify common positions between 3 aligners and store as a list
  Create list of variant objects for each aligner
  For each object in list of variant objects:
    if object's position is in the list of common positions:
      add object to new 'filtered' list

  # Are left with three lists of variant objects, one for each aligner, containing
  # only variants whose positions are found in the common positions list. They
  # are in the same order, as the original variant lists are all ordered lowest
  # to highest position, and that order does not change when filtering the lists.

  Create dictionary of SNPs for the isolate that will contain positions:[ref base, basecall]

  For each equivalent object in the three filtered lists:
    If the DP4 values suggest a REF and agree among all 3 aligners:
      Add position:refbase to SNPs dictionary
    Elif the DP4 values suggest a SNP and agree among all 3 aligners:
      Add position:altbase to SNPs dictionary
    Else: # This is DP4 values 0.20 < x < 0.80, as well as cases where the three
          # aligners do not agree on a call (eg. bwa and novo call SNP, last calls REF)
      Add position:'N' to SNPs dictionary

  # Now we have a dictionary of SNPs and base calls for the isolate, and are ready
  # to create the full sequence:

  Create empty sequence string

  For each high quality position in the common positions master list:
    If position appears in SNPs dictionary for that isolate:
      Append the base call to the sequence
    Elif position does not appear in SNPs dictionary for that isolate:
      Append the reference base to the sequence from the master dictionary of positions/reference bases

  # Now just have to write the sequence to file, as well as create a file
  # containing the SNPs dictionary, which serves as a reference for which positions
  # were called SNP, REF, or N for that isolate.

  If output fasta file exists:
    Open fasta file and append sequence
  Else:
    Create fasta file and append sequence

  Create base calls file for isolate and write SNPs dictionary to it

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Notes about files required for this analysis:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run this script, two sets of vcf files must exist for each isolate: vcf files
with minimum quality 30 to identify common, high quality positions, and vcf files
with minimum quuality 25 for the actual variant calling.

These vcf files should be generated using bcftools filter. The combining of files
from each aligner for each isolate is HARD CODED into this script; it therefore
MUST be run on isolates with EXACTLY THREE (3) sets of vcf files.

This script does not filter the vcf files themselves for quality, depth, etc.,
but does call SNPs based on DP4 values.

Essentially, follow the workflow in my lab notes. Perhaps I shall write out a
detailed workflow for the lab notes at some point.
