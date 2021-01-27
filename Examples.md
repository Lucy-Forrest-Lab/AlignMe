# Examples of how to run AlignMe

AlignMe has a number of different modes. Here we provide examples of various standard alignment types.
1. [Pair-wise Sequence-to-sequence](#sequence\-to\-sequence-alignments)
2. 


## Sequence-to-sequence alignments

This section refers to alignments based on two protein sequences
provided by the user. The sequence properties used for aligning the
sequences can vary. For such an alignment, two sequence files and a
similarity score file must be provided.

The flags for the fasta files are:  
`--fasta_file1 FILE1`  
`--fasta_file2 FILE1`  

The two amino acid sequences that you want to align must be in [fasta
format](Formats.md) and in separate files (FILE1 and FILE2). 

`--similarity_score_file FILE`

After this flag you have to provide filename of a file containing
information about the type of alignment you wish to create. This file
can be set up in 3 different ways depending on the type of alignment:
(a) SequenceSimilarity, (b) ScaleSimilarity or (c) ProfileSimilarity or
(d) any combination of (a), (b) and (c):

5.1.1) Similarity Score File for a pair-wise sequence alignment 
----------------------------------------------------------------

### a) based on similarity matrices:

When using similarity matrices the corresponding line in the
similarity\_score\_file should have the following format:

Weight \<double\> type SequenceSimilarity file \<filename\>

An example of a similarity score file for an alignment based on a
matrix:

![](media/image2.png){width="6.263888888888889in"
height="0.4583333333333333in"}

The word following "type" describes the kind of alignment you wish to
create. In this case SequenceSimilarity creates an alignment based on a
substitution matrix. The filename has to contain the complete address of
the file containing the matrix you want to use (only relative paths are
used in the examples). An example substitution matrix is given in
section 4.1.

### b) Similarity Score File for a pair-wise sequence alignment based on Position Specific Matrices

*Weight \<double\> type PostionSpecificSimilarity PSSM1 \<filename\>
PSSM2 \<filename\>*

An example of a similarity score file for an alignment based on Position
Specific Matrices:

![](media/image3.png){width="6.305555555555555in"
height="0.3194444444444444in"}

The type PositionSpecificSimilarity is used to generated alignments
based on the PSSMs which are provided subsequently. The file of PSSM1
has to be based on the sequence provided by the flag -fasta\_file1 and
the file after PSSM2 has to correspond to the sequence of the flag
--fasta\_file2.

The principle of PositionSpecificSimilarity for the replacement of amino
acid A from sequence 1 by amino acid B from sequence 2 is: From PSSM1
the value for the replacement of A with B is taken, from PSSM2 the value
for the replacement of B with A is taken and then the average of those
values is calculated (their sum divided by 2).

Moreover, there is another method of using PSSMs for alignments
supported by AlignMe called "PositionSpecificSimilarity". The syntax is
similar to those of "PositionSpecificSimilarity":

![](media/image4.png){width="6.305555555555555in"
height="0.3194444444444444in"}

However, the calculation is different. From PSSM1 all 20 values of amino
acid A are compared to the corresponding values of amino acid B in PSSM1
(i.e., likeliness of A-\>A of PSSM1 with B-\> A of PSSM2, A-\>C of 1,
B-\>C of 2 and so on), their differences are summed up and then this
value is divided by 20 (= number of amino acids types).

### c) Similarity Score File for a pair-wise sequence alignment based on scales (e.g. hydrophobicity):

When using scales the corresponding line in the similarity\_score\_file
should have the following format:

*Weight \<double\> type ScaleSimilarity file \<filename\> windowtype
\<string\> windowsize \<integer\>*

An example of a similarity score file for an alignment based on a
matrix:

![](media/image5.png){width="6.263888888888889in"
height="0.3611111111111111in"}

The type ScaleSimilarity is used to create an alignment based on a
scale. The filename refers to the file containing the scale you wish to
use. In such a scale a value is assigned for all 20 amino acid types
(see section 4.2 for format details). The similarity between two amino
acids is calculated as the difference between their scale values.

Currently, 5 different window types are supported in AlignMe. These
windows provide the option to smooth the scale values by averaging the
values over a subsequence. These are: none, rectangular, triangular,
sinoid, and zigzag (see section 5 for detailed information).

### d) Similarity Score File for a pair-wise sequence alignment based on profiles (e.g. secondary structure predictions)

When aligning using profiles, the similarity\_score\_file should have
the following format:

*Weight \<double\> type UniversalProfileSimilarity column \<double\>
headerlines \<double\> profile1 \<filename1\> profile2 \<filename2\>*

An example of a similarity score file for an alignment based on
profiles:

![](media/image6.png){width="6.25in" height="0.20833333333333334in"}

Here, the type UniversalProfileSimilarity is used in order to align
user-specified profiles. A profile contains values in a certain column
corresponding to a certain amino acid of the sequence (see section 6.5
for format). Therefore, the profiles must have the same length as the
sequences. If the length of one of the profiles does not match that of
the corresponding sequence (i.e., profile1 corresponds to fasta\_file1,
and profile2 to fasta\_file2), an error message will be given. In
addition, the column number from which the values will be taken must be
provided. The number given after the tag "headerlines" describes the
number of lines that will be skipped at the beginning of the
profile-file. This option is useful if there are comments or other
information at the beginning of the profile file that you do not want to
include in the alignment.

### e) Similarity Score File for a pair-wise sequence alignment based on a combination of matrices, scales or profiles

AlignMe allows combinations of input types. Each input type is defined
in a separate row of the similarity\_score\_file. The following is an
example of a similarity score file containing a combination of inputs:

![](media/image7.png){width="6.25in" height="0.4027777777777778in"}

5.1.2) Example commands for Pair-wise Alignments
------------------------------------------------

Here, we acquaint you with the basic commands of AlignMe with some
examples using required flags. More information concerning optional
flags is explained in section 4.

Change directory to the main AlignMe folder to run the following
commands.

#### A) Pair-wise alignment of 2 sequences based on a substitution matrix

*./alignme1.2.exe --fasta\_file1 ./examples/fastas/1H2S\_A.fa
--fasta\_file2 ./examples/fastas/2EI4\_A.fa --similarity\_score\_file
./examples/similarity\_scorefiles/matrix.txt*

You will receive the following standard warning messages indicating the
usage of default values:

*No gap opening penalty provided. The default value 10 will be used.*

*No gap extension penalty provided. The default value 1 will be used.*

*No termini extension penalty provided. The default value 1 will be
used.*

*No termini opening penalty provided. The default value will be used. *

*You did not provide a filename for the output of the sequence
alignment. It will be written to aligned\_sequences.aln*

To modify the default values and to define custom output files, see
section 4.

Only these 3 basic input flags that you have used to create this
alignment are required flags. The alignment is now stored in
aligned\_sequences.aln. This alignment has been created based on the
BLOSUM substitution matrix.

#### B) Pair-wise alignment of 2 sequences based on a PositionSpecificSubstitutionMatrix (PSSM)

*/alignme1.2.exe --fasta\_file1 ./examples/fastas/1H2S\_A.fa
--fasta\_file2 ./examples/fastas/2EI4\_A.fa --similarity\_score\_file
./examples/similarity\_scorefiles/PSSM.txt --output\_aligned\_profiles
my\_aligned\_PSSMs.aln *

With this command, an alignment is generated based on
"PositionSpecificSimilarity" and the aligned sequences are stored in the
file "*my\_aligned\_PSSMs.aln*".

However, there is also the scoring type
"ProfilePositionSpecificSimilarity" available, which has been described
in section 3.1 and can be used with the following command:

*/alignme1.2.exe --fasta\_file1 ./examples/fastas/1H2S\_A.fa
--fasta\_file2 ./examples/fastas/2EI4\_A.fa --similarity\_score\_file
./examples/similarity\_scorefiles/PSSMprofile.txt
--output\_aligned\_profiles my\_aligned\_profile PSSMs.aln *

Note that only the file after the flag "*--similarity\_score\_file"* is
different!

#### C) Pair-wise alignment of 2 sequences based on a hydrophobicity scale 

*alignme1.2.exe --fasta\_file1 ./examples/fastas/1H2S\_A.fa
--fasta\_file2 ./examples/fastas/2EI4\_A.fa --similarity\_score\_file
./examples/similarity\_scorefiles/scale.txt*

When creating alignments based on scales or profiles, it can be useful
to use the optional flag *--*output\_aligned\_profiles (see section 4
for detailed information) to create an additional output file containing
the aligned values of each sequence position. For example:

*alignme1.2.exe --fasta\_file1 ./examples/fastas/1H2S\_A.fa
--fasta\_file2 ./examples/fastas/2EI4\_A.fa --similarity\_score\_file
./examples/similarity\_scorefiles/scale.txt *

*--output\_aligned\_profiles my\_aligned\_profiles.aln *

The aligned profiles are now written to my\_aligned\_profiles.aln, while
the sequence alignment is still written to aligned\_sequences.aln. To
get a better overview of the underlying hydrophobicity of your sequence,
the profile file can then be plotted, e.g. with xmgrace or gnuplot.

#### D) Pair-wise alignment of 2 sequences based on secondary-structure predictions

In this example two per-residue predictions are aligned with each other.
These predictions were obtained with the secondary structure predictor
PsiPred, but you can use any kind of program that creates a profile
(i.e., transmembrane predictors, secondary structure predictors etc.).

*alignme1.2.exe --fasta\_file1 ./examples/fastas/1H2S\_A.fa
--fasta\_file2 --similarity\_score\_file
./examples/similarity\_scorfiles/profile.txt --output\_aligned\_profiles
my\_aligned\_profiles.aln *

#### E) Pair-wise alignment of 2 sequences based on combined inputs 

*alignme1.2.exe --fasta\_file1 ./examples/1KPL.fa --fasta\_file2
./examples/1OTS.fa \\\
--similarity\_score\_file ./examples/similarity\_scorefiles/combined.txt
--output\_aligned\_profiles my\_aligned\_profiles.aln *

This alignment is built using three different input types, i.e.
amino-acid substitution (matrix), hydrophobicity (scale) and secondary
structure (profile). For difficult alignments involving sequences with
low sequence similarity, such combinations are usually more accurate
than one input alone.

5.2) Alignment of two Multiple Sequence Alignments 
---------------------------------------------------

This section describes alignments based on two multiple sequence
alignments (MSAs) provided by the user. Currently, MSAs can only be used
with scales. The scale values are averaged to generate an averaged
profile, and then the two profiles are pair-wise aligned. This method is
also referred to as family-averaged profile alignment.

The flags for the input alignments are:

***--msa\_file1 \<filename1\>***

***--msa\_file2 \<filename2\>***

The sequences you want to align must be in fasta format in separate
files (filename1 and filename2). Each of these files must contain one or
more sequences. Each sequence in the file must start with the '\>'
symbol, followed by a header (it can also be left empty). All sequences
of a file need to have the same length (including gaps), because it is
assumed that these sequences are all aligned with each other. At least
one sequence has to be provided per file.

***--fraction\_allowed\_gaps \<double\>***

Each column of the multiple sequence alignment is checked for the
fraction of gaps that it contains. If the fraction of gaps in a given
column is higher than this "fraction of gaps" threshold value, this
column will not be considered in the alignment. Default value = 0.5,
i.e. columns in which more than 50 % of all positions are gaps are
skipped and are not considered in the alignment.

***--similarity\_score\_file \<filename\>***

This flag requires you to provide a file containing information about
the type of alignment you want to do. Currently, only alignments based
on (hydrophobicity) scales, and using a triangular window are supported
for alignment of two averaged multiple sequence alignments.

### 5.2.1) Similarity Score File for an alignment of two MSAs

For aligning two MSAs the line in the similarity\_score\_file should
have the following format:

*Weight \<double\> type ScaleSimilarity file \<filename\> windowtype
msa\_triangular windowsize \<integer\>*

A valid example of a similarity score file:

![](media/image8.png){width="6.25in" height="0.4444444444444444in"}

In the above example the filename refers to the file containing the
scale with which you wish to create an alignment. If amino acids of the
submitted sequences are not in the corresponding scale, an error message
will be given.

In the current version of AlignMe only the triangular\_msa window type
is supported (see section 5 for detailed information about sliding
window types).

The length of the sequence must be longer than the chosen window size;
otherwise the program will quit with an error.

### 5.2.2) Example command

Enter the folder AlignMe main folder and test the following command:

*alignme1.2.exe --msa\_file1 ./examples/bcct.fa --msa\_file2
./examples/deda.fa \\\
--similarity\_score\_file ./examples/simscore\_msa.txt
--fraction\_allowed\_gaps 0.5*

Warnings are shown about default values; these can be ignored. Take a
look at the aligned profiles, which have been written by default to
aligned\_profiles.aln.

5.3) Pair-wise profile-to-profile alignments
--------------------------------------------

This section is about alignments based only on profiles. In contrast to
approaches discussed in sections 2.1 and 2.2, you must not provide an
amino acid sequence, which improves the speed of the alignment.
Moreover, this option allows alignments of any kind of profile and is
therefore not restricted to sequence alignments.

### 5.3.1) Similarity Score File for a pairwise profile-to-profile alignment

*--*similarity\_score\_file \<filename\>

For an alignment without sequences you only can use the type
UniversalProfileSimilarity:

The similarity\_score\_file has to look like:

*Weight \<double\> type UniversalProfileSimilarity column \<double\>
headerlines \<double\> profile1 \<filename1\> profile2 \<filename2\>*

Two profiles have to be provided. A profile contains corresponding
values in a certain column (for more details see section 1) and you have
to choose the column that will be used. Headerlines describes the number
of lines that will be skipped at the beginning of the profile-file. This
option is useful if there are comments or other information at the
beginning of a file you do not want to include for the alignment.

### 5.3.2) Example command

Enter the folder AlignMe main folder and test the following command:

*alignme1.2.exe --similarity\_score\_file
./examples/similarity\_score\_files/profile.txt
--output\_aligned\_profiles my\_aligned\_profiles.aln*

You can have a look at your aligned profiles that will are stored in
my\_aligned\_profiles.aln

