# Running AlignMe on the Web Server

Jump to section on:
- [Getting started](#getting-started)
- [Pairwise Sequence Alignment](#pairwise-sequence-alignment) 
   - Alignment Parameters
   - Batch Mode
   - Anchors
- [Align Multiple Sequence Alignments](#align-multiple-sequence-alignments)   
   - Alignment Parameters
- Examples 
   - Pairwise
   - Batch
   - MSA
   



### Getting started 

To run AlignMe on the command line in a terminal use:  

`alignme`

To get an overview of available options, see [overview](Flag_overview.md), or call:  

`./alignme -help`

The program reads a variety of flags that are explained in
this manual, and that can be provided in any order. Flgas start with the
 '-' symbol and expect a filename or a value to follow. There are
[required flags](#Required-inputs) that have no default values defined and that will cause
error messages when missing. The [optional flags](#Optional-flags) have default values, and
these are printed as warnings when the flag is not given. (Look at the
[examples section](#Examples.md) for more information about warnings.)

---

### Pairwise Sequence Alignment

There are three required inputs for every AlignMe calculation:
1) the first [sequence or profile](#sequence-or-profile-inputs)
2) the second [sequence or profile](#sequence-or-profile-inputs)
3) a file containing the similarity metrics to be used, or ["similarity score file"](#similarity-score-file)

These inputs can be combined to run AlignMe in a number of different modes. General information about the required inputs is given below. Specific examples of various standard alignment types are provided here:
1. [Pair-wise Sequence-to-sequence](#sequence\-to\-sequence-alignments). These are the most common.
2. [Alignment of two Multiple Sequence Alignments](#Alignment-of-two-Multiple-Sequence-Alignments)
3. [Pair-wise profile-to-profile alignments](#Pair\-wise-profile\-to\-profile-alignments)

---

### Align Multiple Sequence Alignments

