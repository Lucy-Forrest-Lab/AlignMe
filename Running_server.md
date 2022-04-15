# Running AlignMe on the Web Server

Jump to section on:
- [Getting started](#getting-started)
- [Pairwise Sequence Alignment](#pairwise-sequence-alignment) 
   - [Alignment Inputs and Parameters](#alignment-inputs-and-parameters)
   - [Anchors](#anchors)
   - [Batch Mode](#batch-mode)
- [Align Multiple Sequence Alignments](#align-multiple-sequence-alignments)   
   - Alignment Parameters
- Examples 
   - Pairwise
   - Batch
   - MSA
   



### Getting started 

To access the AlignMe Web Server go to: 

https://www.bioinfo.mpg.de/AlignMe/



---

### Pairwise Sequence Alignment

#### Alignment Inputs and Parameters
There are three required inputs for every pairwise AlignMe calculation:
1) the first sequence (can be uploaded or pasted into the specified box)
2) the second sequence (can be uploaded or pasted in to the specified box)
3) selection of an optimized parameter mode (or user-defined parameters)

These inputs can be combined to run AlignMe in a number of different modes. There are four modes with optimized parameters (including gap opening and extension parameters) that the user can choose:
1. Fast: based on a substituion matrix (VTML) and hydrophobicity scale (HWvH)
2. P Mode: based on a Position Specific Substitution Matrix (PSSM) which is generated during the AlignMe run 
3. PS Mode: based on a PSSM and secondary structure prediction which are both generated during the AlignMe run
4. PST Mode: based on a PSSM, secondary structure prediction, and a transmembrane prediction all of which are generated during the AlignMe run

In addition to these four modes, a user can also choose to define their own parameters. This includes choosing gap penalties and choosing a PSSM, secondary structure prediction, and/or a transmembrane prediction to be generated during the AlignMe run. The user can also upload their own PSSM, secondary structure predictions, and/or transmembrane prediction. 

#### Anchors 

AlignMe has the option to include position restraints, known as anchors, in the sequence alignment. The user needs to specify which amino acids to be matched in each sequence as well as an anchor constraint strength factor, known as an anchor weight.

There are two ways an anchor can be input on the AlignMe Web Server:
1. Input on the server page
2. Uploaded in a file

If the user chooses to upload a file containing the anchor positions, each anchor has to be given as a line in the anchor file, following the scheme described in the [Formats section](#Formats.md):  

`position_in_first_sequence   position_in_second_sequence   strength`

e.g.: 
`25  36 1000`

This would align residue 25 in the first sequence with residue 36 of the second sequence with a relative strength of 1000. 

More detail can be found in the section on [Anchors](#Anchors.md).

#### Batch Mode

The AlignMe Web Server also allows the user to run many pairwise alignments at once, this is called Batch Mode. In order to use Batch Mode, instead of pasting individual sequences into each input box, the user can paste all of the sequences to be aligned into the box, or upload a file with multiple FASTA sequences. For the provided sets of sequences, every sequence in the first input will be aligned with every sequence in the second input, up to a maximum of 1000 total alignments per submission. Every alignment will be pairwise and will use the same set of parameters. The result page will provide the alignment and plots (based on which mode is selected) for the first set of sequences that were aligned. The rest of the results can be downloaded via links on the results page. 

Batch mode also allows the input of anchors. To use anchors in Batch Mode, the user must upload a file (cannot use the input boxes on the server) in the following format: 

`position of sequence in input file: position_in_first_sequence   position of sequence in input file: position_in_second_sequence   strength`

e.g.:
`1:18 7:48 1000`

This would align residue 18 in the first sequence with residue 48 in the seventh sequence with a relative strength of 1000.


---

### Align Multiple Sequence Alignments

