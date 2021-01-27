# Examples of how to run AlignMe

AlignMe has a number of different modes. Here we provide examples of various standard alignment types.
1. [Pair-wise Sequence-to-sequence](#sequence\-to\-sequence-alignments)
2. [Alignment of two Multiple Sequence Alignments](#Alignment-of-two-Multiple-Sequence-Alignments)
3. [Pair-wise profile-to-profile alignments](#Pair\-wise-profile\-to\-profile-alignments)


## Sequence-to-sequence alignments

This section refers to alignments based on two protein sequences
provided by the user. The sequence properties used for aligning the
sequences can vary. For such an alignment, two sequence files and a
similarity score file must be provided.

The flags for the fasta files are:  
`-fasta_file1 FILE1`  
`-fasta_file2 FILE1`  

The two amino acid sequences that you want to align must be in [fasta
format](Formats.md) and in separate files (FILE1 and FILE2). 

`-similarity_score_file FILE`

After this flag you have to provide filename of a file containing
information about the type of alignment you wish to create. This file
can be set up in 3 different ways depending on the type of alignment:
1. [Sequence Similarity Matrices](#Using-similarity-matrices), 
2. [Scale Similarity](#Using-Scale-similarity) or 
3. [Profile Similarity](#Using-Profile-similarity) or 
4. any [combination](#Combinations-of-inputs) of (a), (b) and (c).

##### Using Similarity Matrices
When using similarity (aka substitution) matrices the corresponding line in the
similarity_score_file should have the following format:  

`weight: <double> type: SequenceSimilarity file: <filename>`

for example:

`weight: 1.0 type: SequenceSimilarity file: ./examples/matrices/blosum62.dat`

The word following "type" describes the kind of alignment you wish to
create. In this case **SequenceSimilarity** creates an alignment based on a
substitution matrix. The filename has to contain the complete address of
the file containing the matrix you want to use (only relative paths are
used in the examples). Example substitution matrices are available [here](https://github.com/Lucy-Forrest-Lab/AlignMe/tree/master/examples/matrices/).

b) Similarity Score File for a pair-wise sequence alignment based on Position Specific Matrices

`weight: <double> type: PostionSpecificSimilarity PSSM1: <filename> PSSM2: <filename>`

An example of a similarity score file for an alignment based on Position Specific Matrices:

`weight; 1.0 type: PostionSpecificSimilarity PSSM1: ./examples/PSSMs/1H2S.pssm PSSM2: ./examples/PSSMs/2EI4.pssm`

The type **PositionSpecificSimilarity** is used to generated alignments
based on the PSSMs which are provided subsequently. The file of PSSM1
has to be based on the sequence provided by the flag -fasta_file1 and
the file after PSSM2 has to correspond to the sequence of the flag
-fasta_file2.

The principle of PositionSpecificSimilarity for the replacement of amino
acid A from sequence 1 by amino acid B from sequence 2 is: From PSSM1
the value for the replacement of A with B is taken, from PSSM2 the value
for the replacement of B with A is taken and then the average of those
values is calculated (their sum divided by 2).

Moreover, there is another method of using PSSMs for alignments
supported by AlignMe called "PositionSpecificSimilarity". The syntax is
similar to those of "PositionSpecificSimilarity":
`weight: 1.0 type: ProfilePostionSpecificSimilarity PSSM1: ./examples_best/1KPL.pssm PSSM2: ./examples_best/1OTS.pssm`

However, the calculation is different. From PSSM1 all 20 values of amino
acid A are compared to the corresponding values of amino acid B in PSSM1
(i.e., likeliness of A->A of PSSM1 with B->A of PSSM2, A->C of 1,
B->C of 2 and so on), their differences are summed up and then this
value is divided by 20 (= number of amino acid types).

### Using Scale Similarity
Similarity Score File for a pair-wise sequence alignment based on scales (e.g. hydrophobicity):

When using scales the corresponding line in the similarity_score_file
should have the following format:

`weight: <double> type: ScaleSimilarity file: <filename> windowtype: <string> windowsize: <integer>`

An example of a similarity score file for an alignment based on a matrix:

`weight: 1.0 type: ScaleSimilarity file: ./examples/scales/KD.txt windowtype: triangular windowsize: 13`

The type ScaleSimilarity is used to create an alignment based on a
scale. The filename refers to the file containing the scale you wish to
use. In such a scale a value is assigned for all 20 amino acid types
(see [section on format details]()). The similarity between two amino
acids is calculated as the difference between their scale values.

Currently, 5 different window types are supported in AlignMe. These
windows provide the option to smooth the scale values by averaging the
values over a subsequence. These are: none, rectangular, triangular,
sinoid, and zigzag (see section 5 for detailed information).

#### Using Profile Similarity
Similarity Score File for a pair-wise sequence alignment based on profiles (e.g. secondary structure predictions)

When aligning using profiles, the similarity score file should have the following format:

`weight: <double> type: UniversalProfileSimilarity column: <double> headerlines: <double> profile1: <filename1> profile2: <filename2>`

An example of a similarity score file for an alignment based on profiles:  

`weight: 1.0 type: UniversalProfileSimilarity column: 5 headerlines: 1 profile1: ./examples/profiles/1H2S_A.ss2 profile2: ./examples/profiles/2EI4_A.ss2`

Here, the type UniversalProfileSimilarity is used in order to align
user-specified profiles. A profile contains values in a certain column
corresponding to a certain amino acid of the sequence (see [section on formats]()). 
Therefore, the profiles must have the same length as the
sequences. If the length of one of the profiles does not match that of
the corresponding sequence (i.e., profile1 corresponds to fasta_file1,
and profile2 to fasta_file2), an error message will be given. In
addition, the column number from which the values will be taken must be
provided. The number given after the tag "headerlines" describes the
number of lines that will be skipped at the beginning of the
profile-file. This option is useful if there are comments or other
information at the beginning of the profile file that you do not want to
include in the alignment.

### Combinations of Inputs
Similarity Score File for a pair-wise sequence alignment based on a combination of matrices, scales or profiles

AlignMe allows combinations of input types. Each input type is defined
in a separate row of the similarity_score_file. The following is an
example of a similarity score file containing a combination of inputs:

```
weight: 1.0 type: SequenceSimilarity file: ./examples/matrices/blosum62.mat
weight: 1.0 type: ScaleSimilarity file: ./examples/scales/KD.txt windowtype: triangular windowsize: 13
weight: 1.0 type: UniversalProfileSimilarity column: 5 headerlines: 1 profile1: ./examples/profiles/1H2S_A.ss2 profile2: ./examples/profiles/2EI4_A.ss2
```

### Example commands for Pair-wise Alignments

Here, we acquaint you with the basic commands of AlignMe with some
examples using required flags. More information concerning optional
flags is explained [here]().

Change directory to the main AlignMe folder to run the following commands.

#### Pair-wise alignment of 2 sequences based on a substitution matrix

```./alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa -similarity_score_file ./examples/similarity_scorefiles/matrix.txt```

You will receive the following standard warning messages indicating the
usage of default values:

- No gap opening penalty provided. The default value 10 will be used.
- No gap extension penalty provided. The default value 1 will be used.
- No termini extension penalty provided. The default value 1 will be
used.
- No termini opening penalty provided. The default value will be used. 
- You did not provide a filename for the output of the sequence
alignment. It will be written to **aligned_sequences.aln**

To modify the default values and to define custom output files, see
[sections on inputs](#Running.md) and [overview of flags](#Flag_overview.md).

Only these 3 basic input flags that you have used to create this
alignment are required flags. The alignment is now stored in
**aligned_sequences.aln**. This alignment has been created based on the
BLOSUM substitution matrix.

#### Pair-wise alignment of 2 sequences based on a PositionSpecificSubstitutionMatrix (PSSM)

```
./alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa 
-similarity_score_file ./examples/similarity_scorefiles/PSSM.txt -output_aligned_profiles my_aligned_PSSMs.aln
```

With this command, an alignment is generated based on
"PositionSpecificSimilarity" and the aligned sequences are stored in the
file **my_aligned_PSSMs.aln**.

However, there is also the scoring type
"ProfilePositionSpecificSimilarity" available, which has been described
in [section](3.1) and can be used with the following command:

```
./alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa 
-similarity_score_file ./examples/similarity_scorefiles/PSSMprofile.txt -output_aligned_profiles my_aligned_profile_PSSMs.aln
```

Note that only the file after the flag `-similarity_score_file` differs from the previous example.

#### Pair-wise alignment of 2 sequences based on a hydrophobicity scale 

```
alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa 
-similarity_score_file ./examples/similarity_scorefiles/scale.txt
```

When creating alignments based on scales or profiles, it can be useful
to use the optional flag --output_aligned_profiles (see sections on [flags]() and [outputs]()
for detailed information) to create an additional output file containing
the aligned values of each sequence position. For example:

```
alignme.exe -fasta\_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa 
-similarity_score_file ./examples/similarity_scorefiles/scale.txt -output_aligned_profiles my_aligned_profiles.aln
```

The aligned profiles are now written to **my_aligned_profiles.aln**, while
the sequence alignment is still written to **aligned_sequences.aln**. To
get a better overview of the underlying hydrophobicity of your sequence,
the profile file can then be plotted, e.g. with xmgrace or gnuplot.

#### Pair-wise alignment of 2 sequences based on secondary-structure predictions

In this example two per-residue predictions are aligned with each other.
These predictions were obtained with the secondary structure predictor
PsiPred, but you can use any kind of program that creates a profile
(i.e., transmembrane predictors, secondary structure predictors etc.).

```
alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa 
-similarity_score_file  ./examples/similarity_scorefiles/profile.txt -output_aligned_profiles my_aligned_profiles.aln
```

####Pair-wise alignment of 2 sequences based on combined inputs 

```
alignme.exe -fasta_file1 ./examples/1KPL.fa -fasta_file2 ./examples/1OTS.fa 
-similarity_score_file ./examples/similarity_scorefiles/combined.txt -output_aligned_profiles my_aligned_profiles.aln
```

This alignment is built using three different input types, i.e.
amino-acid substitution (matrix), hydrophobicity (scale) and secondary
structure (profile). For difficult alignments involving sequences with
low sequence similarity, such combinations are usually more accurate
than one input alone.

## Anchoring a pairwise alignment
To introduce fixed positions in the alignment, you can call a list of anchors using the following syntax:
```
alignme -similarity_score_file <similarityscorefile> -fasta_file1 <filename> -fasta_file2 <filename> 
-anchors <anchorsfile>
```

---

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

## Pair-wise profile-to-profile alignments

This section is about alignments based only on profiles. In contrast to
approaches discussed in sections 2.1 and 2.2, you must not provide an
amino acid sequence, which improves the speed of the alignment.
Moreover, this option allows alignments of any kind of profile and is
therefore not restricted to sequence alignments.

### Similarity Score File for a pairwise profile-to-profile alignment

`-similarity_score_file <filename>`

For an alignment without sequences you only can use the type
UniversalProfileSimilarity:

The similarity_score_file has to look like:

`weight: <double> type: UniversalProfileSimilarity column: <double> headerlines: <double> profile1: <filename1> profile2: <filename2>`

Two profiles have to be provided. A profile contains corresponding
values in a certain column (for more details see [section on something]()) and you have
to choose the column that will be used. Headerlines describes the number
of lines that will be skipped at the beginning of the profile-file. This
option is useful if there are comments or other information at the
beginning of a file you do not want to include for the alignment.

##### Example command

Enter the folder AlignMe main folder and test the following command:

```
alignme.exe -similarity_score_file ./examples/similarity_score_files/profile.txt 
-output_aligned_profiles  my_aligned_profiles.aln
```

You should have a look at your aligned profiles that will be written to **my_aligned_profiles.aln**.



