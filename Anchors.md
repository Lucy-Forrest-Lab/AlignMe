# Anchoring positions in AlignMe alignments

Alignments can be optimized by including expert knowledge. Any number of anchors can be defined for positions in pairs of sequences that are known to belong together. 

Each anchor has to be given as a line in the anchor file, following the scheme described in the [Formats section](#Formats.md):  
`position_in_first_sequence   position_in_second_sequence   strength`

e.g.: 
`25  36 1000`

This would align residue 25 in the first sequence with residue 36 of the second sequence with a relative strength of 1000. The strength is subject to preference and may need some try and error adjustment. However, we generally recommend the use of a weight that universally imposes the matching of the two residues; specifically, weight values > 10 will typically be sufficient to match the two positions, no matter the complexity or sequence similarity of the sequences. To remove all uncertainty, the default value given to the weight is 1000. However, by gradually increasing the weight in the range of 0.1 to 100 and assessing the impact on the alignment (i.e. whether the two positions are matched, and what is the effect on the local region on the alignment), the optimum value can be identified for a given pair of sequences.
 
The anchor file is then called using the following flag:

`-anchors \<filename\>`
