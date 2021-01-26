# Introduction to AlignMe

AlignMe is a program that allows the user to perform pair-wise alignments for
- two sequences 
- two multiple sequence alignments
- two profiles

The similarity (*Sim(i,j)*) between two amino acids (*i* and *j*) may be
measured by their evolutionary relationship in terms of a substitution
matrix (*S(a,b)*) or by differences of biochemical profiles, such as
hydrophobicity or predicted structural values, computed as (*\|V(a) -- V(b)\|*).
AlignMe can combine any number of input profiles (*n*)
with any number of input substitution matrices (*m*) and can use different weights
(*w*) to guide the combination of those inputs in generating the alignment. 

This similarity measure is calculated for all amino acid pairs of two
sequences in the Needleman-Wunsch algorithm. Additionally the user can
choose between several gap penalty schemes for the treatment of
mismatching stretches (section 4.2).

In all modes of usage, the user has to provide two files, containing the
sequences, and a so-called similarity-score file, containing the
user-defined composition of the similarity scores (sections 3 and 5 ).

The user can also define any number of constraints or anchors. The anchors can be tuned 
in strength from preference to enforcement.

