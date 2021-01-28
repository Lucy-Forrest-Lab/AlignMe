# Running AlignMe

Jump to section on:
- [Getting Started](#getting-started)
- [Required Inputs](#required-inputs)
   - [Sequence Inputs](#sequence-inputs)
   - [Similarity Score File](#similarity-score-file)
- [Optional Flags](#optional-flags)  
   - [Setting Anchors](#anchors)
   - [Defining gap penalties](#gap-penalties)  
   - [Alignment algorithm](#alignment-algorithm)  
- [Output Files](Output.md)
- [Examples](Examples.md)
   - [Pair-wise Sequence-to-sequence](#sequence\-to\-sequence-alignments). 
   - [Alignment of two Multiple Sequence Alignments](#Alignment-of-two-Multiple-Sequence-Alignments)
   - [Pair-wise profile-to-profile alignments](#Pair\-wise-profile\-to\-profile-alignments)
- [File Formats](Formats.md)
- [Overview of all flags](Flag_overview.md)


### Getting started 

To run AlignMe on the command line in a terminal use:  

`alignme.exe`

To get an overview of available options, see [overview](Flag_overview.md), or call:  

`./alignme.exe -help`

The program reads a variety of flags that are explained in
this manual, and that can be provided in any order. Flgas start with the
 '-' symbol and expect a filename or a value to follow. There are
[required flags](#Required-inputs) that have no default values defined and that will cause
error messages when missing. The [optional flags](#Optional-flags) have default values, and
these are printed as warnings when the flag is not given. (Look at the
[examples section](#Examples.md) for more information about warnings.)

---

### Required inputs

There are three required inputs for every AlignMe calculation:
1) the first [sequence or profile](#sequence-or-profile-inputs)
2) the second [sequence or profile](#sequence-or-profile-inputs)
3) a file containing the similarity metrics to be used, or ["similarity score file"](#similarity-score-file)

These inputs can be combined to run AlignMe in a number of different modes. General information about the required inputs is given below. Specific examples of various standard alignment types are provided here:
1. [Pair-wise Sequence-to-sequence](#sequence\-to\-sequence-alignments). These are the most common.
2. [Alignment of two Multiple Sequence Alignments](#Alignment-of-two-Multiple-Sequence-Alignments)
3. [Pair-wise profile-to-profile alignments](#Pair\-wise-profile\-to\-profile-alignments)


#### Sequence or profile inputs
To parse the sequence inputs, one of the following two pairs of flags has to be provided:  

`-fasta_file1 <filename1> -fasta_file2 <filename2>`  
or  
`-msa_file1 <filename1> -msa_file2 <filename2>`  
fasta_file = filename of the file containing a primary sequence in Fasta format.
msa_file = filename of multiple-sequence alignment that will be converted to a family-averaged hydropathy profile.

#### Similarity score file
The so-called "similarity score file" contains information about the similarity metrics that you wish to
use to construct your alignment. The file is called from the command line using the following flag:

`-similarity_score_file <filename>` 

The similarity score metrics can be any combination of primary sequence, hydropathy scale, or profiles such as from secondary structure or transmembrane predictions. Each row in the similarity score file must start with a weight, which
indicates the influence of the chosen parameter on the alignment. If you are using only one type of similarity measure, then this value should be 1.0. 
The entry after the option "type:" describes the type of alignment to be created.
In the simplest example, for  an alignment based on a substitution matrix, the similarity score file would contain this line:

`weight: 1.0 type: SequenceSimilarity file: ./examples/matrices/blosum62.mat`

Every input line must end with a carriage return/new line, otherwise the entry will not be read. 

More details about similarity score files can be found in the [examples](#Examples.md).

---

### Optional flags

#### Anchors
Adding a list of positions that will be fixed can be carried out using this flag:
`-anchors <filename>`
The format of this file is described [here](#Formats.md). More detail can be found in the section on [Anchors](#Anchors.md).

#### Gap penalties 
In any alignment two conditions exist, namely, sections of sequences
match or they do not match. Mismatches are reflected by the introduction
of gaps. Therefore, the quality of the alignment will be influenced not
only by the similarity measures but also the criteria of introducing
gaps; the latter are known as gap penalties. There are no "correct"
values for gap penalties, instead they have to be adjusted by experience
or optimization based on knowledge of "correct" alignments.
Correspondingly, AlignMe has a variety of flags to control gap
penalties.

The basic gap penalty scheme consists of two different penalties:  

`-gap_opening_penalty <value>`  
`-gap_extension_penalty <value>`  

Typically, alignments with fewer and longer gaps are preferred over many
short gaps. Therefore to open a new gap is usually assigned a higher
penalty than to extend an existing gap. By default a gap opening penalty
of 10 and a gap extension penalty of 1 are used. Different values can be
set using these flags. Both values must be positive integers or
fractions (no commas). This basic scheme corresponds to the "2 penalties" option on the webserver.

This basic scheme can be extended by the use of 2 additional gap
penalties, which control gaps at the ends (N- and C-termini) of the
sequences:  

`-termini_gap_opening_penalty <value>`  
`-termini_gap_extension_penalty <value>` 

The latter flags allow gaps at the end of a sequence to be treated
differently from gaps within the sequence. If the lengths of the
sequences differ significantly it can be useful to assign lower
penalties to termini regions. If these values are not set, then by
default they are assigned the same values as `gap_open_penalty` and
`gap_extension_penalty`. Again, both values must be positive integers or
fractions (no commas). This scheme corresponds to the "4 penalties" option on the webserver, and is used in the AlignMe P mode.

Finally, one may also wish to penalize gap penalties in the core of conserved regions, such as membrane-spanning segments or secondary-structure elements. To this end, a more advanced penalty scheme can be introduced using the flags:  

`-below_threshold_gap_opening_penalty <value>`  
`-below_threshold_gap_extension_penalty <value>`  
`-above_threshold_gap_opening_penalty <value>`  
`-above_threshold_gap_extension_penalty <value>`  
`-thresholds_for_penalties <value> <value> ...`  

With the last flag one can define thresholds for the scales or profiles, above or below which different penalties are applied.
The other four flags provide the penalties for those different regions. For example: if a
secondary structure prediction determines a helix probability above 0.5
one can avoid the introduction of gaps in that helix by defining a
threshold of 0.5 and by providing larger values for the 'above'
penalties than for the 'below' penalties.  

Note that the values passed to the threshold flag have to match the
definitions in the similarity score file. For example, the flag:  

`-thresholds_for_penalties 0.5 0.5`  

will apply these thresholds and the according gap penalties to the first two scales/profiles defined
in the similarity score file, e.g. helix and hydrophobicity
probabilities. Any scale/profile/matrix defined later in the similarity
score file, (e.g. coil probability and hydrophobicity) will not be
subject to this penalty scheme. In addition, if these five flags are
provided, both termini gap penalties must also be provided. 
This scheme corresponds to the "6 penalties" option on the webserver, and is used in AlignMe Fast, PS, and PST modes.

#### Alignment Algorithm 
The user has the option to change the algorithm being used for the
alignment. Currently, only one option is available, i.e.
"global_affine", which is also the default value.

`-algorithm <name>`  
