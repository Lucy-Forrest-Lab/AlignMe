# Anchoring positions in AlignMe alignments

Alignments can be optimized by including expert knowledge. Any number of anchors can be defined for positions in pairs of sequences that are known to belong together. 

Each anchor has to be given as a line in the anchor file, following the scheme described in the [Formats section](#Formats.md):  
`position_in_first_sequence   position_in_second_sequence   strength`

e.g.: 
`25  36 1000`

This would align residue 25 in the first sequence with residue 36 of the second sequence with a relative strength of 1000. The strength is subject to preference and may need some try and error adjustment. A strength of 1000 should enforce the anchor for most cases, since the magnitude of the input data in scales, matrices and profiles is typically <10. 
 
The anchor file is then called using the following flag:

`-anchors \<filename\>`
