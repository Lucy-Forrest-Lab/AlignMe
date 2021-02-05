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
//!
//!
//!
//! @author: Rene Staritzbichler
//! @date: 9.6.2011
/////////////////////////////////////////////////////////////////////////

#ifndef ALIGNMENT_VARIABLES_H_
#define ALIGNMENT_VARIABLES_H_

#include "function.t.h"
#include <set>

struct AlignNut
{
	std::vector< std::pair< int, int> >    m_Alignment;
	int                                    m_FirstID;
	int 								   m_SecondID;
	double								   m_LastElement;
};

struct AlignmentVariables
{
	  double
		  gap_opening_penalty,
		  gap_extension_penalty,
		  termini_gap_opening_penalty,
		  termini_gap_extension_penalty,
		  below_threshold_gap_opening_penalty,
		  above_threshold_gap_opening_penalty,
		  below_threshold_gap_extension_penalty,
		  above_threshold_gap_extension_penalty,
//		  normalization_factor,
		  last_element,                              // ??
		  fraction_allowed_gaps;

	  size_t
		  line_length;  // for the output of the alignments

	  std::string
		  base,                                        // ??
		  first_fasta_id,
		  second_fasta_id,
		  algorithm,
		  first_fasta_file,
		  second_fasta_file,
		  out_file,
		  similarity_score_file,
		  outputfile_aligned_sequences,
		  outputfile_aligned_ids,
		  outputfile_aligned_profiles,
		  first_msa_file,
		  second_msa_file,
		  alignment_output_format,
		  gap_value;

	  ShPtr< Function< std::vector< double>, double> >
		  gap_opening_penalty_function,
		  gap_extension_penalty_function;


	  std::vector< double>
		  thresholds;

	  std::vector< std::string>
		  file_endings,
		  write_profiles_header1,
		  write_profiles_header2;

	  // A GeneralizedAminoAcid has as initializations-list m_Type and m_Profiles
	  std::vector< Sequence>
		  first_sequences,                    //  ??
		  second_sequences,                   //  ??
		  first_msa,
		  second_msa;

	  // creates structure for the scores which will be a pair of AA combined with a certain value and the number
	  // of the corresponding shared pointer...
	  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
		  scores,
		  alignment_scores;

	  // creates structure for alignment with a number combined with two other numbers... ???
	  std::multimap< double, AlignNut>
		  alignments;

	  std::set<char>
		  defined_amino_acids;

	  AlignNut
		  alignment_data;

	  std::pair< double, std::vector< std::pair< int, int> > >
		  score_and_alignment;


	  Sequence
		  first_sequence,
		  second_sequence,
		  first_averaged_alignment,
		  second_averaged_alignment;

	  // for averaged MSA
	  std::vector< int>
		  first_ignored_positions,
		  second_ignored_positions;

	  Matrix< DynamicProgrammingMatrixElement>
		  dynamic_programing_matrix;

	  std::vector< double>
		  weights;

	  std::vector< Triplet<int,int,double> >
	  	  anchors;
};

void SetDefaults( AlignmentVariables &NUT)
{
    NUT.gap_opening_penalty = 10.0;
    NUT.gap_extension_penalty = 1.0;
//    NUT.normalization_factor = 1.0;
	NUT.algorithm = "global_affine";
	NUT.out_file = "AlignMe1.1_";
	NUT.line_length = 60;
	NUT.outputfile_aligned_sequences = "aligned_sequences.aln";
	NUT.outputfile_aligned_ids = "aligned_ids.aln";
	NUT.outputfile_aligned_profiles = "";
	NUT.first_msa_file = "";
	NUT.second_msa_file = "";
	NUT.alignment_output_format =  "clustalw";
	NUT.fraction_allowed_gaps = 0.5;
	NUT.gap_value = "?0";
}


#endif /* ALIGNMENT_VARIABLES_H_ */
