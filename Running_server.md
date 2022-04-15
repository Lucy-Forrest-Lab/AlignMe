# Running AlignMe on the Web Server

Jump to section on:
- [Getting started](#getting-started)
- [Pairwise Sequence Alignment](#pairwise-sequence-alignment) 
   - [Alignment Inputs and Parameters](#alignment-inputs-and-parameters)
   - Anchors
   - Batch Mode
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
1) the first sequence
2) the second sequence
3) selection of an optimized parameter mode (or user-defined parameters)

These inputs can be combined to run AlignMe in a number of different modes. There are four modes with optimized parameters (including gap opening and extension parameters) that the user can choose:
1. Fast: based on a substituion matrix (VTML) and hydrophobicity scale (HWvH)
2. P Mode: based on a Position Specific Substitution Matrix (PSSM) which is generated during the AlignMe run 
3. PS Mode: based on a PSSM and secondary structure prediction which are both generated during the AlignMe run
4. PST Mode: based on a PSSM, secondary structure prediction, and a transmembrane prediction all of which are generated during the AlignMe run

In addition to these four modes, a user can also choose to define their own parameters. This includes choosing gap penalties and choosing a PSSM, secondary structure prediction, and/or a transmembrane prediction to be generated during the AlignMe run. The user can also upload their own PSSM, secondary structure predictions, and/or transmembrane prediction. 

#### Anchors 


---

### Align Multiple Sequence Alignments

