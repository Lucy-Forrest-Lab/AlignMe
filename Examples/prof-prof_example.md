## Pair-wise profile-to-profile alignments

This section is about alignments based only on profiles. In contrast to
approaches discussed in sections 2.1 and 2.2, you must not provide an
amino acid sequence, which improves the speed of the alignment.
Moreover, this option allows alignments of any kind of profile and is
therefore not restricted to sequence alignments.

### Similarity Score File for a pairwise profile-to-profile alignment

`-similarity_score_file <filename>`

For an alignment without sequences you only can use the type
UniversalProfileSimilarity:

The similarity_score_file has to look like:

`weight: <double> type: UniversalProfileSimilarity column: <double> headerlines: <double> profile1: <filename1> profile2: <filename2>`

Two profiles have to be provided. A profile contains corresponding
values in a certain column (for more details see [section on something]()) and you have
to choose the column that will be used. Headerlines describes the number
of lines that will be skipped at the beginning of the profile-file. This
option is useful if there are comments or other information at the
beginning of a file you do not want to include for the alignment.

##### Example command

Enter the folder AlignMe main folder and test the following command:

```
alignme.exe -similarity_score_file ./examples/similarity_score_files/profile.txt 
-output_aligned_profiles  my_aligned_profiles.aln
```

You should have a look at your aligned profiles that will be written to **my_aligned_profiles.aln**.


