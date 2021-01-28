## Alignment of two multiple sequence alignments 

This section describes alignments based on two multiple sequence
alignments (MSAs) provided by the user. Currently, MSAs can only be used
with hydrophobicity scales. The scale values at each position in the each MSA are averaged to generate two "family-averaged" profiles, and then the two profiles are pair-wise aligned. This method is
also referred to as family-averaged profile alignment.

In addition to the [required flags](./Running.md#required-inputs), which includes:  
- the filenames of the multiple-sequence alignments in [fasta format](Formats.md), with all sequences having the same length (including gaps), and  
- the [similarity score file](#Similarity-score-file-for-an-alignment-of-two-MSAs),   

the following flag is required:

`-fraction_allowed_gaps <double>`

Each column of the multiple sequence alignment is checked for the
fraction of gaps that it contains. If the fraction of gaps in a given
column is higher than this "fraction of gaps" threshold value (in the range of 0 to 1), this
column will not be considered in the alignment. Default value = 0.5,
i.e. columns in which more than 50 % of all positions are gaps are
skipped and are not considered in the alignment.

The command for launching an example MSA-MSA alignment is provided [below](#Example-command)

---

### Similarity score file for an alignment of two MSAs

For aligning two MSAs the line in the similarity score file should have the following format:  
`weight: <double> type: ScaleSimilarity file: <filename> windowtype: msa_triangular windowsize: <integer>`

A valid example of a similarity score file:

![](media/image8.png)

In the above example the filename refers to the file containing the
scale with which you wish to create an alignment. If amino acids of the
submitted sequences are not in the corresponding scale, an error message
will be given.

In the current version of AlignMe only the triangular_msa window type
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
