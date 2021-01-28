## Sequence-to-sequence alignments

This section refers to alignments based on two protein sequences
provided by the user. The sequence properties used for aligning the
sequences can vary. For such an alignment, two sequence files and a
similarity score file must be provided.

The flags for the fasta files are:  
`-fasta_file1 <filename1>`  
`-fasta_file2 <filename2>`  

The two amino acid sequences that you want to align must be in [fasta
format](Formats.md) and in separate files (FILE1 and FILE2). 

`-similarity_score_file <filename>`

After this flag you have to provide filename of a file containing information about the type of alignment you wish to create. This file can be set up in 3 different ways depending on the type of alignment. The similarity score file syntax for each of these options can be found here:
1. [Sequence similarity matrices](#Using-similarity-matrices)
2. [Position-specific similarity matrices](#Using-position\-specific-matrices) or
3. [Scales for similarity](#Using-scales-for-similarity) or 
4. [Profile similarity](#Using-profiles-for-similarity) or 
5. any [combination](#Combinations-of-inputs) of (a), (b) and (c).
The linked examples below illustrate the file inputs required to run each of these types of alignment, starting with the syntax for the similarity score file for each of the above 5 options. 

The commands to [run](#launching-pairwise-alignments) these examples can be found here:
1. [Sequence similarity matrices](#Using-similarity-matrices)
2. [Position-specific similarity matrices](#Using-position\-specific-matrices) or
3. [Scales for similarity](#Using-scales-for-similarity) or 
4. [Profile similarity](#Using-profiles-for-similarity) or 
5. any [combination](#Combinations-of-inputs) of (a), (b) and (c).


We also provide pre-optimized gap penalties and input files for four different modes: [P](), [PS](), [PST](), and [Fast]() modes.

---

### Using similarity matrices
When using similarity (aka substitution) matrices the corresponding line in the
similarity score file should have the following format:  

`weight: <double> type: SequenceSimilarity file: <filename>`

for example:  

`weight: 1.0 type: SequenceSimilarity file: ./examples/matrices/blosum62.dat`

The word following "type" describes the kind of alignment you wish to create. In this case **SequenceSimilarity** creates an alignment based on a
substitution matrix. Example substitution matrices are available [here](https://github.com/Lucy-Forrest-Lab/AlignMe/tree/master/examples/matrices/).

---

### Using position-specific matrices
A similarity score file for a pair-wise sequence alignment based on position-specific substitution matrices (PSSMs) requires the following format:

`weight: <double> type: PostionSpecificSimilarity PSSM1: <filename> PSSM2: <filename>`

for example:

`weight: 1.0 type: PostionSpecificSimilarity PSSM1: ./examples/PSSMs/1H2S.pssm PSSM2: ./examples/PSSMs/2EI4.pssm`

The type **PositionSpecificSimilarity** is used to generated alignments based on PSSMs. The file of PSSM1 has to be based on the sequence provided by the flag -fasta_file1 and the file after PSSM2 has to correspond to the sequence of the flag -fasta_file2.

The principle of PositionSpecificSimilarity for the replacement of amino acid A from sequence 1 by amino acid B from sequence 2 is: From PSSM1 the value for the replacement of A with B is taken, from PSSM2 the value for the replacement of B with A is taken and then the average of those values is calculated (their sum divided by 2). 

Moreover, there is another method of using PSSMs for alignments supported by AlignMe called **ProfilePositionSpecificSimilarity**. The syntax is
similar to those of "PositionSpecificSimilarity":  
`weight: 1.0 type: ProfilePostionSpecificSimilarity PSSM1: ./examples_best/1KPL.pssm PSSM2: ./examples_best/1OTS.pssm`

However, the calculation is different. From PSSM1 all 20 values of amino acid A are compared to the corresponding values of amino acid B in PSSM1 (i.e., likeliness of A->A of PSSM1 with B->A of PSSM2, A->C of 1, B->C of 2 and so on), their differences are summed up and then this value is divided by 20 (= number of amino acid types). 

---

### Using scales for similarity
Similarity Score File for a pair-wise sequence alignment based on scales (e.g. hydrophobicity):

When using scales the corresponding line in the similarity_score_file should have the following format:

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

---

### Using profiles for similarity
When aligning using profiles (e.g. secondary structure predictions), the similarity score file should have the following format:

`weight: <double> type: UniversalProfileSimilarity column: <double> headerlines: <double> profile1: <filename1> profile2: <filename2>`

An example of a similarity score file for an alignment based on profiles:  

```weight: 1.0 type: UniversalProfileSimilarity column: 5 headerlines: 1 profile1: ./examples/profiles/1H2S_A.ss2 profile2: ./examples/profiles/2EI4_A.ss2```

Here, the type UniversalProfileSimilarity is used in order to align user-specified profiles. A profile contains values in a certain column corresponding to a certain amino acid of the sequence (see [section on formats](Formats.md)). Therefore, the profiles must have the same length as the sequences. If the length of one of the profiles does not match that of the corresponding sequence (i.e., profile1 corresponds to fasta_file1, and profile2 to fasta_file2), an error message will be given. In addition, the column number from which the values will be taken must be provided. The number given after the tag "headerlines" describes the number of lines that will be skipped at the beginning of the profile-file. This option is useful if there are comments or other information at the beginning of the profile file that you do not want to include in the alignment.

---

### Combinations of inputs
Creating a pair-wise sequence alignment based on a combination of matrices, scales or profiles is a powerful feature of AlignMe. The standard Fast, PS, and PST modes use combinations of inputs. Each input type is defined in a separate row of the similarity_score_file. The following is an example of a similarity score file containing a combination of inputs:

```
weight: 1.0 type: SequenceSimilarity file: ./examples/matrices/blosum62.mat
weight: 1.0 type: ScaleSimilarity file: ./examples/scales/KD.txt windowtype: triangular windowsize: 13
weight: 1.0 type: UniversalProfileSimilarity column: 5 headerlines: 1 profile1: ./examples/profiles/1H2S_A.ss2 profile2: ./examples/profiles/2EI4_A.ss2
```

---
---

## Launching pairwise alignments
The above example provides the syntax for creating the similarity score file. The next section explains how to call those input files when running AlignMe.  More information concerning optional flags is explained in the [Flag overview](#Flag_overview.md) and [Running AlignMe](#Running.md) sections.

The following examples can be run by first changing directory to the main AlignMe folder, or by copying the Examples/ folder and its contents to your working directory.
The examples are:

1. [Pairwise alignment of 2 sequences based on a substitution matrix](#Pairwise-alignment-of-2-sequences-based-on-a-substitution-matrix)
2. [Pairwise alignment of 2 sequences based on a position specific substitution matrix](#Pairwise-alignment-of2-sequences-based-on-a-position-specifis-substitution-matrix)
3. [Pairwise alignment of 2 sequences based on a hydrophobicity scale](#Pairwise-alignment-of-2-sequences-based-on-a-hydrophobicity-scale)
4. [Pairwise alignment of 2 sequences based on secondary-structure predictions](#Pairwise-alignment-of-2-sequences-based-on-secondary\-structure-predictions)
5. [Pairwise alignment of 2 sequences based on combined inputs](#Pairwise-alignment-of-2-sequences-based-on-combined-inputs)
6. [Anchoring a pairwise alignment](#Anchoring-a-pairwise-alignment)


#### Pairwise alignment of 2 sequences based on a substitution matrix
The simplest alignment uses the BLOSUM substitution matrix to align two sequences

```
./alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa \
              -similarity_score_file ./examples/similarity_scorefiles/matrix.txt
```

You will receive the following standard warning messages indicating the
usage of default values:  
- No gap opening penalty provided. The default value 10 will be used.
- No gap extension penalty provided. The default value 1 will be used.
- No termini extension penalty provided. The default value 1 will be used.
- No termini opening penalty provided. The default value will be used. 
- You did not provide a filename for the output of the sequence alignment. It will be written to **aligned_sequences.aln**

To modify the default values and to define custom output files, see [sections on inputs](#Running.md) and [overview of flags](#Flag_overview.md).
This example illustrates usage of the three required input flags, and creates an alignment in the file called **aligned_sequences.aln**.

#### Pairwise alignment of 2 sequences based on a position specific substitution matrix
At the next level of complexity, we compare two sequences using a PSSM.
```
./alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa \
              -similarity_score_file ./examples/similarity_scorefiles/PSSM.txt \
              -output_aligned_profiles my_aligned_PSSMs.aln
```

With this command, an alignment is generated based on "PositionSpecificSimilarity" and the aligned sequences are stored in the
file **my_aligned_PSSMs.aln**. However, a second scoring type, called "ProfilePositionSpecificSimilarity" is also available, which has been described
in [section](3.1) and can be used with the following command:

```
./alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa \
              -similarity_score_file ./examples/similarity_scorefiles/PSSMprofile.txt \
              -output_aligned_profiles my_aligned_profile_PSSMs.aln
```

Note that only the file specified by the flag `-similarity_score_file` differs from the previous example.

#### Pairwise alignment of 2 sequences based on a hydrophobicity scale 
For aligning two sequences based on a hydrophobicity scale, the basic command is as follows:
```
./alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa \
              -similarity_score_file ./examples/similarity_scorefiles/scale.txt
```

When creating alignments based on scales or profiles, it can be useful
to use the optional flag --output_aligned_profiles (see sections on [flags](#Flag_overview.md) and [outputs](#Outputs.md)
for detailed information) to create an additional output file containing
the aligned values of each sequence position. For example:

```
alignme.exe -fasta\_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa 
        -similarity_score_file ./examples/similarity_scorefiles/scale.txt \
        **-output_aligned_profiles my_aligned_profiles.aln**
```

The aligned profiles are now written to **my_aligned_profiles.aln**, while
the sequence alignment is still written to **aligned_sequences.aln**. To
get a better overview of the underlying hydrophobicity of your sequence,
the profile file can then be plotted, e.g. with xmgrace or gnuplot.

#### Pairwise alignment of 2 sequences based on secondary-structure predictions
For aligning two per-residue predictions with each other.
These predictions were obtained with the secondary structure predictor
PsiPred, but you can use any kind of program that creates a profile
(i.e., transmembrane predictors, secondary structure predictors etc.).

```
alignme.exe -fasta_file1 ./examples/fastas/1H2S_A.fa -fasta_file2 ./examples/fastas/2EI4_A.fa \
            -similarity_score_file  ./examples/similarity_scorefiles/profile.txt \
            -output_aligned_profiles my_aligned_profiles.aln
```

#### Pairwise alignment of 2 sequences based on combined inputs 

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
