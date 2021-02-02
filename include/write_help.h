/////////////////////////////////////////////////////////////////////////
//!
//  This file is part of the AlignMe program
//
//  Copyright (C) 2010 by Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//  AlignMe@rzg.mpg.de
//
//  AlignMe is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  AlignMe is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//!
//!
//!
//! @author: Rene Staritzbichler, Marcus Stamm, Kamil Khafizov
//! @date: 9.6.2011
/////////////////////////////////////////////////////////////////////////


#ifndef WRITE_HELP_H_
#define WRITE_HELP_H_


void WriteHelpForAlignPairs()
{
	  std::cout << "\n";
	  std::cout << "This is the usage information for AlignMe Version 1.2\n";
	  std::cout << "\n";
	  std::cout << "Flags for input files:\n";
	  std::cout << "-similarity_score_file [file]         file containing information about the type \n";
	  std::cout << "                                      of alignment you want to do. See the\n";
	  std::cout << "                                      manual for detailed information!\n";
	  std::cout << "-fasta_file1 [file]                   file containing an amino acid sequence\n";
	  std::cout << "                                      in fasta format \n";
	  std::cout << "-fasta_file2 [file]                   file containing an amino acid sequence\n";
	  std::cout << "                                      in fasta format \n";
	  std::cout << "-msa_file1 [file]                     file containing a multiple sequence\n";
	  std::cout << "                                      in which all sequence are the\n";
      std::cout << "                                      same length\n";
	  std::cout << "-msa_file2  [file]                    file containing a multiple sequence\n";
	  std::cout << "                                      in which all sequence are the\n";
      std::cout << "                                      same length\n";
	  std::cout << "\n";
	  std::cout << "-similarity_score_file must always be provided!\n";
	  std::cout << "-fasta_file1 and -fasta_file2 must be provided together \n";
	  std::cout << " The same is true for -msa_file1 and msa_file2.\n";
	  std::cout << "\n";
	  std::cout << "Flags for output files:\n";
	  std::cout << "-output_aligned_sequences [file]      file to which the aligned amino acid\n";
	  std::cout << "                                      sequences are printed\n";
	  std::cout << "-output_aligned_profiles [file]       file to which the aligned profile values\n";
	  std::cout << "                                      are printed \n";
	  std::cout << "-extract_from_MSA_sequences_with_ids  extract a sequence from the averaged\n";
	  std::cout << " [value][value]                       multiple sequence alignment. Value1 is\n";
	  std::cout << "                                      the i-th sequence of msa_file1 and\n";
	  std::cout << "                                      Value2 is the j-th sequence of msa_file2\n";
	  std::cout << "-output_extracted_sequences [file]    file to which the sequences that are \n";
	  std::cout << "                                      extracted from the averaged MSA are \n";
	  std::cout << "                                      written to\n";
	  std::cout << "\n";
	  std::cout << "You can choose one of the following gap penalty sets\n";
	  std::cout << "A combination is not possible!\n";
	  std::cout << "\n";
	  std::cout << "Flags for a set of 2 gap penalties:\n";
	  std::cout << "-gap_opening_penalty [value]          penalty for opening gaps\n";
	  std::cout << "-gap_extension_penalty [value]        penalty for extending gaps \n";
	  std::cout << "\n";
	  std::cout << "Flags for a set of 4 gap penalties:\n";
	  std::cout << "-gap_opening_penalty [value]          penalty for opening gaps within\n";
	  std::cout << "                                      the sequence \n";
	  std::cout << "-gap_extension_penalty [value]        penalty for extending gaps within\n";
	  std::cout << "                                      the sequence \n";
	  std::cout << "-termini_gap_opening_penalty          additional penalty for opening gaps\n";
	  std::cout << " [value]                              at the end of the sequence\n";
	  std::cout << "-termini_gap_extension_penalty        additional penalty for extending gaps\n";
	  std::cout << " [value]                              at the end of the sequence\n";
	  std::cout << "\n";
	  std::cout << "Flags for a set of 6 gap penalties:	\n";
	  std::cout << "-thresholds_for_penalties [values]     defines a threshold according to which\n";
	  std::cout << "                                       gap penalties are assigned to a given\n";
	  std::cout << "                                       position of the sequence\n";
	  std::cout << "-below_threshold_gap_opening_penalty   penalty for opening gaps opposite to \n";
	  std::cout << " [value]                               residues with a value below the\n";
	  std::cout << "                                       chosen threshold\n";
	  std::cout << "-below_threshold_gap_extension_penalty penalty for extending gaps opposite to \n";
	  std::cout << " [value]                               residues with a value below the\n";
	  std::cout << "                                       chosen threshold\n";
	  std::cout << "-above_threshold_gap_opening_penalty   penalty for opening gaps opposite to \n";
	  std::cout << " [value]                               residues with a value above the\n";
	  std::cout << "                                       chosen threshold\n";
	  std::cout << "-above_threshold_gap_extension_penalty penalty for extending gaps opposite to \n";
	  std::cout << " [value]                               residues with a value above the\n";
	  std::cout << "                                       chosen threshold\n";
	  std::cout << "-termini_gap_opening_penalty           additional penalty for opening gaps\n";
	  std::cout << " [value]                               at the end of the sequence\n";
	  std::cout << "-termini_gap_extension_penalty         additional penalty for extending gaps\n";
	  std::cout << " [value]                               at the end of the sequence\n";
	  std::cout << "\n";
	  std::cout << "Other flags:";
	  std::cout << "\n";
	  std::cout << "-fraction_allowed_gaps [value]         columns of the MSA in which more than \n ";
	  std::cout << "                                       [value] positions are gaps are skipped  \n";
	  std::cout << "                                       and not considered in the alignment\n";
	  std::cout << "-profile_gap_value_for_plotting        define a value to be assigned to gaps\n";
	  std::cout << " [value]                               in the output aligned profiles \n";
	  std::cout << "-anchors [filename]                    file containing positions in both sequences\n";
	  std::cout << "                                       to be aligned, each line in the file is one anchor:\n";
	  std::cout << "                                       the first line contains the number of anchors used:\n";
	  std::cout << "                                       e.g.:  \n";
	  std::cout << "                                       1              (nr of anchors)\n";
	  std::cout << "                                       25  36   2.5   (pos in seq A, pos in seq B, strength)\n";
	  std::cout << "\n";
	  std::cout << "\n";
	  exit (-1);
}


#endif /* WRITE_HELP_H_ */
