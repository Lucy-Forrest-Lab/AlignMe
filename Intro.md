1) Introduction
===============

AlignMe is a program that allows the user to perform pair-wise
alignments for

-   two sequences (section 3.1),

-   two multiple sequence alignments (section 3.2) two profiles (section
    3.3)

Similarity (*Sim~i,j~*) between two amino acids (*I* and *J*) may be
measured by their evolutionary relationship in terms of a substitution
matrix (*S~A,B~*) or by differences of biochemical profiles, such as
hydrophobicity or predicted structural values

(*\|V~a~ -- V~b~\|*). AlignMe can combine any number of profiles (*n*)
with any number of substitution matrices (*m*) with different weights
(*w*) to guide the alignment (section 2.1.1)

This similarity measure is calculated for all amino acid pairs of two
sequences in the Needleman-Wunsch algorithm. Additionally the user can
choose between several gap penalty schemes for the treatment of
mismatching stretches (section 4.2).

In all modes of usage, the user has to provide two files, containing the
sequences, and a so-called similarity-score file, containing the
user-defined composition of the similarity scores (sections 3 and 5 ).

New in Version 1.2:

-   The user can define any number of constraints or anchors. The
    anchors can be tuned in strength from preference to enforcement
    (section 4.4).

-   Pairwise alignments can be performed in a batch mode, which is most
    relevant for the server

To get an overview of available options call:

`*./alignme.exe -help*`


