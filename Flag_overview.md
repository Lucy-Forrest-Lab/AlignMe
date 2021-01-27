## Overview of all flags

### Flags for input files:

  Flag                                   | Description
  -------------------------------------- | ------------------------------------------------------------------------------------------
  --similarity_score_file FILE           | file containing information about the type of alignment to be run
  --fasta_file1 FILE                     | file containing an amino acid sequence in fasta format
  --fasta_file2 FILE                     | file containing an amino acid sequence in fasta format
  --msa_file1 FILE                       | file containing a multiple sequence alignment, with entries all the same length
  --msa_file2 FILE                       | file containing a multiple sequence alignment, with entries all the same length

Note the following requirements:  
  `--similarity_score_file` must always be provided.  
  
  `--fasta_file1` and `--fasta_file2` must be provided together.  
  The same is true for `--msa_file1` and `--msa_file2`.


### Flags for output files:

Flag                               | Description
-----------------------------------| -----------------------------------
 --output_aligned_sequences FILE   | file to which the aligned amino acid sequences are printed
 --output_aligned_profiles         | file to which the aligned profile values are printed     
 --extract_from_MSA_sequences_with_ids  \[value1\] \[value2\]    | from each averaged multiple sequence alignment, a sequence is extracted. Value1 = i-th sequence from msa_file1. Value2 = j-th sequence from msa_file2
 --output_extracted_sequences FILE  | file to which the sequences that are extracted from the averaged MSA are written


## Gap penalty settings

You can choose one of the following gap penalty sets. A combination is not possible.

### Flags for a set of 2 gap penalties:

  Flag                                            | Description
  ----------------------------------------------- | ----------------------------
  --gap_opening_penalty \[value\]                 | penalty for opening gaps
  --gap_extension_penalty \[value\]               | penalty for extending gaps


### Flags for a set of 4 gap penalties:

  Flag                               | Description
  ----------------------------------------------- | -----------------------------------------------------------------
  --gap_opening_penalty \[value\]                 |  penalty for opening gaps within the sequence
  --gap_extension_penalty \[value\]               |  penalty for extending gaps within the sequence
  --termini_gap_opening_penalty \[value\]         |  additional penalty for opening gaps at the end of the sequence
  --termini_gap_extension_penalty \[value\]       |  additional penalty for extending gaps at the end of the sequence

### Flags for a set of 6 gap penalties:

  Flag                               | Description
  ------------------------------------------------- | -------------------------------------------------------------------------------------------------------
  --thresholds_for_penalties \[values\]             |  defines a threshold according to which gap penalties are assigned
  --below_threshold_gap_opening_penalty \[value\]   |  penalty for opening gaps opposite residues with a value below the  threshold
  --below_threshold_gap_extension_penalty \[value\] |  penalty for extending gaps opposite residues with a value below the  threshold
  --above_threshold_gap_opening_penalty  \[value\]  |  penalty for opening gaps opposite residues with a value above the  threshold
  --above_threshold_gap_extension_penalty \[value\] |  penalty for extending gaps opposite residues with a value above the  threshold
  --termini_gap_opening_penalty \[value\]           |  additional penalty for opening gaps at the end of the sequence
  --termini_gap_extension_penalty \[value\]         |  additional penalty for extending gaps at the end of the sequence


### Other flags

  Flag                               | Description
  -------------------------- | -----------
  --anchors FILE             | filename of the file containing a list of anchors, e.g. anchors.txt
  --fraction_allowed_gaps \[value\]           | columns of the MSA in which more than \[value\] positions are gaps are skipped and not considered in the alignment
  --profile_gap_value_for_plotting \[value\]  | define a value to be assigned to gaps in the output aligned profiles; default is ?0
  --alignment_output_format                   | Formatting type of the output alignment. Allowed types are ClustalW and fasta. By default, alignments are written in ClustalW format.
