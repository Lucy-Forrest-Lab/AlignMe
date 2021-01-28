## Alignment of two Multiple Sequence Alignments 

This section describes alignments based on two multiple sequence
alignments (MSAs) provided by the user. Currently, MSAs can only be used
with scales. The scale values are averaged to generate an averaged
profile, and then the two profiles are pair-wise aligned. This method is
also referred to as family-averaged profile alignment.

The flags for the input alignments are:

`-msa_file1 <filename1>`  
`-msa_file2 <filename2>`

The sequences you want to align must be in fasta format in separate
files (filename1 and filename2). Each of these files must contain one or
more sequences. Each sequence in the file must start with the '\>'
symbol, followed by a header (it can also be left empty). All sequences
of a file need to have the same length (including gaps), because it is
assumed that these sequences are all aligned with each other. At least
one sequence has to be provided per file.

`-fraction_allowed_gaps <double>`

Each column of the multiple sequence alignment is checked for the
fraction of gaps that it contains. If the fraction of gaps in a given
column is higher than this "fraction of gaps" threshold value, this
column will not be considered in the alignment. Default value = 0.5,
i.e. columns in which more than 50 % of all positions are gaps are
skipped and are not considered in the alignment.

`-similarity_score_file <filename>`

This flag requires you to provide a file containing information about
the type of alignment you want to do. Currently, only alignments based
on (hydrophobicity) scales, and using a triangular window are supported
for alignment of two averaged multiple sequence alignments.

### Similarity Score File for an alignment of two MSAs

For aligning two MSAs the line in the similarity score file should
have the following format:  
`Weight <double> type ScaleSimilarity file <filename> windowtype msa_triangular windowsize <integer>`

A valid example of a similarity score file:

![](media/image8.png)

In the above example the filename refers to the file containing the
scale with which you wish to create an alignment. If amino acids of the
submitted sequences are not in the corresponding scale, an error message
will be given.

In the current version of AlignMe only the triangular\_msa window type
is supported (see section 5 for detailed information about sliding
window types).

The length of the sequence must be longer than the chosen window size;
otherwise the program will quit with an error.

### Example command

Enter the folder AlignMe main folder and test the following command:

```
alignme.exe -msa_file1 ./examples/bcct.fa -msa_file2 ./examples/deda.fa 
-similarity_score_file ./examples/simscore_msa.txt -fraction_allowed_gaps 0.5
```

Warnings are shown about default values; these can be ignored. Take a
look at the aligned profiles, which have been written by default to
**aligned_profiles.aln**.
