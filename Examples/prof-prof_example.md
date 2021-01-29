## Pairwise profile-to-profile alignments

This section provides examples of alignments based only on 1-dimensional profiles. In contrast to the pairwise sequence-to-sequence or pairwise MSA-to-MSA alignments, an amino acid sequence is not required, which improves the speed of the alignment. Moreover, this option allows alignments of any kind of profile and is therefore not restricted to sequence alignments.

The first section below describes the expected content of the [similarity score file](#Similarity-score-file-for-a-pairwise-profile\-to\-profile-alignment).  
The command for launching an example profile-profile alignment is provided [below](#Example-command).

---

### Similarity score file for a pairwise profile-to-profile alignment

`-similarity_score_file <filename>`

For an alignment without sequences you only can use the type
UniversalProfileSimilarity:

The similarity_score_file has to look like:

`weight: <double> type: UniversalProfileSimilarity column: <double> headerlines: <double> profile1: <filename1> profile2: <filename2>`

Two profiles have to be provided. A profile contains corresponding
values in a certain column (for more details see [section on profile formats](https://github.com/Lucy-Forrest-Lab/AlignMe/blob/gh-pages/Formats.md#Profiles)) and you must choose the column that will be used. "Headerlines" describes the number
of lines that will be skipped at the beginning of the profile file. This
option is useful if there are comments or other information at the
beginning of a file that you do not want to include for the alignment.

---

##### Example command

Change directory to the main AlignMe folder, or copy the Examples/ folder to your working directory and test the following command:

```
alignme.exe -similarity_score_file ./examples/similarity_score_files/profile.txt \
            -output_aligned_profiles  my_aligned_profiles.aln
```

You should have a look at your aligned profiles, which will be written to **my_aligned_profiles.aln**.
