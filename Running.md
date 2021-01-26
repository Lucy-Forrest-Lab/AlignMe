# Running AlignMe

To get an overview of available options call:
`./alignme.exe --help`

3. [Required Flags]
4. [Output Files](Output.md)
5. [Using Anchors](Anchors.md)
6. [Examples](Examples.md)
7. [File Formats](Formats.md)
8. [Overview of all flags](Flag_overview.md)


3) Required flags 
------------------

Start the program by calling the alignme1.2.exe executable from a
terminal. The program offers a variety of flags that are explained in
this manual, and that can be provided in any order. They start with the
usual '--' symbol and expect a filename or a value to follow. There are
required flags that have no default values defined and that will cause
error messages when missing. Optional flags have default values, and
these are printed as warnings when the flag is not given. Look at
section 5.1.2 to get more information about warnings. 

All modes of AlignMe require a so-called "similarity score file". This
file contains information about the similarity metrics that you wish to
use to construct your alignment. These can be any combination of
sequence, scale, or profile similarity.

Each row in the similarity score file must start with a weight, which
describes the influence of the chosen parameter on the alignment. If you
are using only one type of similarity measure, then this value should be
1.0. An example of a similarity score file for an alignment based on a
substitution matrix is:
`weight: 1.0 type: SequenceSimilarity file: ./examples/matrices/blosum62.mat`

More details about similarity score files are provided [below].

Together with the \'*-similarity\_score\_file FILE*\' flag,one of the following
two pairs of flags has to be provided:

`--fasta\_file1 FILE1 --fasta\_file2 FILE2`

or

`--msa\_file1 FILE1 --msa\_file2 FILE2`

Each of these files should contain one or more amino acid sequences in FASTA format.

4) Optional flags
=================

4.1) Alignment Algorithm 
-------------------------

**--*algorithm \<name\>***

The user has option to change the algorithm being used for the
alignment. Currently, only one option is possible, i.e.
"global\_affine", which is also the default value.

4.2) Gap penalty schemes 
-------------------------

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

**--*gap\_opening\_penalty \<value\>***

**--*gap\_extension\_penalty \<value\>***

Typically, alignments with fewer and longer gaps are preferred over many
short gaps. Therefore to open a new gap is usually assigned a higher
penalty than to extend an existing gap. By default a gap opening penalty
of 10 and a gap extension penalty of 1 are used. Different values can be
set using these flags. Both values must be positive integers or
fractions (no commas).

This basic scheme can be extended by the use of 2 additional gap
penalties, which control gaps at the ends (N- and C-termini) of the
sequences:

**--*termini\_gap\_opening\_penalty \<value\>***

**--*termini\_gap\_extension\_penalty \<value\>***

The latter flags allow gaps at the end of a sequence to be treated
differently from gaps within the sequence. If the lengths of the
sequences differ significantly it is probably useful to assign lower
penalties to termini regions. If these values are not set, then by
default they are assigned the same values as gap\_open\_penalty and
gap\_extension\_penalty. Again, both values must be positive integers or
fractions (no commas).

Finally, an advanced penalty scheme can be introduced using the flags:

**--*below\_threshold\_gap\_opening\_penalty \<value\>***

**--*below\_threshold\_gap\_extension\_penalty \<value\>***

**--*above\_threshold\_gap\_opening\_penalty \<value\>***

**--*above\_threshold\_gap\_extension\_penalty \<value\>***

**--*thresholds\_for\_penalties \<value\> \<value\> \.....***

With the last flag one can define thresholds for the scales or profiles.
The other four flags then allow gaps to be distinguished if they are to
be opened in regions where the profile values are above or below the
thresholds. This scheme can be used to prevent gap insertion in
secondary structure elements or transmembrane regions. For example: if a
secondary structure prediction determines a helix probability above 0.5
one can avoid the introduction of gaps in that helix by defining a
threshold of 0.5 and by providing larger values for the 'above'
penalties than for the 'below' penalties.

Note that the values passed to the threshold flag have to match the
definitions in the similarity score file. For example:
'**--**thresholds\_for\_penalties 0.5 0.5' will apply these thresholds
and the according gap penalties to the first two scales/profiles defined
in the similarity score file, e.g. helix and hydrophobicity
probabilities. Any scale/profile/matrix defined later in the similarity
score file, (e.g. coil probability and hydrophobicity) will not be
subject to this penalty scheme. In addition, if these five flags are
provided, both termini gap penalties must also be provided.

