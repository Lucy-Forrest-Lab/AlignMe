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
//! Collection of functions simplifying alignment creation.
//!
//!
//!
//! @author: Rene Staritzbichler
//! @date: 9.6.2011
/////////////////////////////////////////////////////////////////////////


#ifndef ALIGNMENT_FUNCTIONS_H_
#define ALIGNMENT_FUNCTIONS_H_


#include "needleman_wunsch_affine_gaps.h"


void
AlignPairsGlobalAffine( AlignmentVariables &VARS)
{
	  NeedlemanWunschAffineGaps
		needleman_wunsch_affine_gaps
		(
				VARS.gap_opening_penalty_function,
				VARS.gap_extension_penalty_function,
				VARS.termini_gap_opening_penalty,
				VARS.termini_gap_extension_penalty,
				VARS.first_sequence,
				VARS.second_sequence,
				VARS.scores,
				VARS.dynamic_programing_matrix
		);

//	  std::cout << __FUNCTION__ << " calc matrix" << std::endl;
	  DebugWrite( "needleman wunsch affine gaps: calculate matrix:...");
	  needleman_wunsch_affine_gaps.CalculateMatrix();

//	  std::cout << __FUNCTION__ << " trace back " << std::endl;
	  DebugWrite( "needleman wunsch affine gaps: trace back:...");
	  VARS.score_and_alignment = needleman_wunsch_affine_gaps.TraceBack();

//	  const Matrix< DynamicProgrammingMatrixElement> * matrix_ptr = &needleman_wunsch_affine_gaps.GetMatrix();
//	  std::vector< double> values = ( *matrix_ptr)( matrix_ptr->GetNumberOfRows() - 1, matrix_ptr->GetNumberOfColumns() - 1).GetAffinePathWays();
//	  last_element = *std::max_element( values.begin(), values.end());
//	  DebugWrite ("alignment score is " <<  VARS.score_and_alignment.first);
}


void
AlignPairsGlobal( AlignmentVariables &VARS)
{

	  NeedlemanWunsch
		needleman_wunsch
		(
				VARS.gap_opening_penalty,
				VARS.gap_extension_penalty,
				VARS.termini_gap_opening_penalty,
				VARS.termini_gap_extension_penalty,
				VARS.first_sequence,
				VARS.second_sequence,
				VARS.scores,
				VARS.dynamic_programing_matrix
		);

	  DebugWrite( "needleman wunsch: calculate matrix:...");
	  needleman_wunsch.CalculateMatrix();

	// needleman_wunsch.Write( write);

	  DebugWrite( "needleman wunsch: traceback...");
	  VARS.score_and_alignment = needleman_wunsch.TraceBack();
}

void
AlignPairsLocal( AlignmentVariables &VARS)
{
	SmithWaterman
	smith_waterman
	(
			VARS.gap_opening_penalty,
			VARS.gap_extension_penalty,
			VARS.first_sequence,
			VARS.second_sequence,
			VARS.scores,
			VARS.dynamic_programing_matrix
	);

  DebugWrite( "smith waterman: calculate matrix:...");
  smith_waterman.CalculateMatrix();

  DebugWrite( "smith waterman: traceback...");
  VARS.score_and_alignment = smith_waterman.TraceBack();
}



#endif /* ALIGNMENT_FUNCTIONS_H_ */

