8) Overview of all flags
------------------------

**Flags for input files:**

  -------------------------------------------------- ------------------------------------------------------------------------------------------
  ***--similarity\_score\_file \[file\]        ***   file containing information about the type of alignment you want to do
  ***--fasta\_file1 \[file\]***                      file containing an amino acid sequence in fasta format
  ***--fasta\_file2 \[file\]***                      file containing an amino acid sequence in fasta format
  ***--msa\_file1 \[file\]***                        file containing a multiple sequence alignment in which all sequences are the same length
  ***--msa\_file2 \[file\]***                        file containing a multiple sequence alignment in which all sequences are the same length
  -------------------------------------------------- ------------------------------------------------------------------------------------------

--similarity\_score\_file must always be provided

--fasta\_file1 and --fasta\_file2 must be provided together

The same is true for --msa\_file1 and msa\_file2.\
\
**Flags for output files:**

+-----------------------------------+-----------------------------------+
| ***--output\_aligned\_sequences   | file to which the aligned amino   |
| \[file\] ***                      | acid sequences are printed        |
+-----------------------------------+-----------------------------------+
| ***--output\_aligned\_profiles    | -   file to which the aligned     |
| \[file\] ***                      |     profile values are printed    |
+-----------------------------------+-----------------------------------+
| ***-extract\_from\_MSA\_sequences | -   from each averaged multiple   |
| \_with\_ids                       |     sequence alignment a sequence |
| \[value1\] \[value2\]***          |     is extraced. Value1 is i-th   |
|                                   |     sequence from msa\_file1 and  |
|                                   |     value2 the j-th sequence from |
|                                   |     msa\_file2                    |
+-----------------------------------+-----------------------------------+
| ***-output\_extracted\_sequences  | file to which the sequences that  |
| \[file\]***                       | are extracted from the averaged   |
|                                   | MSA are written to                |
+-----------------------------------+-----------------------------------+

You can choose one of the following gap penalty sets.

A combination is not possible!

**Flags for a set of 2 gap penalties:**

  ----------------------------------------------- ----------------------------
  ***--gap\_opening\_penalty \[value\]      ***   penalty for opening gaps
  ***--gap\_extension\_penalty \[value\] ***      penalty for extending gaps
  ----------------------------------------------- ----------------------------

**Flags for a set of 4 gap penalties:**

  -------------------------------------------------------- ------------------------------------------------------------------
  ***--gap\_opening\_penalty \[value\]      ***            penalty for opening gaps within the sequence
  ***--gap\_extension\_penalty \[value\] ***               penalty for extending gaps within the sequence
  ***--termini\_gap\_opening\_penalty \[value\]      ***   additional penalty for opening gaps at the end of the sequence
  ***--termini\_gap\_extension\_penalty \[value\] ***      additional penalty for extending gaps at the end of the sequence
  -------------------------------------------------------- ------------------------------------------------------------------

**Flags for a set of 6 gap penalties:**

  -------------------------------------------------------------- -------------------------------------------------------------------------------------------------------
  ***--thresholds\_for\_penalties \[values\]***                  defines a threshold according to which gap penalties are assigned to a given position of the sequence
  ***--below\_threshold\_gap\_opening\_penalty \[value\] ***     penalty for opening gaps opposite to residues with a value below the chosen threshold
  ***--below\_threshold\_gap\_extension\_penalty \[value\]***    penalty for extending gaps opposite to residues with a value below the chosen threshold
  ***--above\_threshold\_gap\_opening\_penalty  \[value\] ***    penalty for opening gaps opposite to residues with a value above the chosen threshold
  ***--above\_threshold\_gap\_extension\_penalty \[value\] ***   penalty for extending gaps opposite to residues with a value above the chosen threshold
  ***--termini\_gap\_opening\_penalty \[value\]***               additional penalty for opening gaps at the end of the sequence
  ***--termini\_gap\_extension\_penalty \[value\] ***            additional penalty for extending gaps at the end of the sequence
  -------------------------------------------------------------- -------------------------------------------------------------------------------------------------------

**Other flags**

  ------------------------------------------------------ -------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ***--fraction\_allowed\_gaps \[value\]***              columns of the MSA in which more than \[value\] positions are gaps are skipped and not considered in the alignment
  ***--profile\_gap\_value\_for\_plotting \[value\]***   define a value to be assigned to gaps in the output aligned profiles, the standard value for gaps is ?0
  ***- alignment\_output\_format***                      Formatting type of the alignment to be generated. Allowed types are: ClustalW and fasta. If this flag is not provided, alignments are written in ClustalW format.
  ------------------------------------------------------ -------------------------------------------------------------------------------------------------------------------------------------------------------------------
