## Output files

To specify the filename for an output pairwise sequence alignment, use the following flag:

`--output_aligned_sequences <filename>`

By default the alignment is written to a file called **aligned_sequences.aln**.


To specify the filename for an output pairwise profile alignment (also known as a family-averaged hydropathy alignment or msa-msa alignment), use the following flag:

`--output_aligned_profiles <filename>`

By default the alignment will be written to a file called **aligned_profiles.aln**.


To specify the filename for any aligned profiles generated during the alignment, use the following flag:

`--profile_gap_value_for_plotting <value>`

By default gaps will be written as a "?0" in the output profile alignment. This option is particularly useful for plotting your aligned
sequences including gaps. The resultant file can be plotted e.g., with gnuplot or, after reformating with the attached perl script `reformat_profile_output.pl`, with xmgrace/qtgrace.

