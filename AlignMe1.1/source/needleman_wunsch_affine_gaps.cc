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
//!  Needleman Wunsch algorithm with affine gap penalties
//! (different penalties for opening and continuing gaps).
//! As described in "Biological sequence analysis"
//! from Durbin, Eddy, Krogh, Mitchison
//!
//!
//! @author: Marcus Stamm, Kamil Khafizov,
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#include "../include/needleman_wunsch_affine_gaps.h"
#include <limits>
#include <climits>


// default constructor
NeedlemanWunschAffineGaps::NeedlemanWunschAffineGaps()
  : m_Matrix( 0, 0)
{}


// construct from data
// needleman_wunsch( gap_opening_penalty, gap_extension_penalty, termini_gap_opening_penalty, termini_gap_extension_penalty,
// first_sequence, second_sequence, scores);
// assign these values directly to class members with an initialization-list

NeedlemanWunschAffineGaps::NeedlemanWunschAffineGaps
( 
		 const ShPtr< Function< std::vector< double>, double> > &GAP_OPENING_PENALTY_FCT,
		 const ShPtr< Function< std::vector< double>, double> > &GAP_EXTENSION_PENALTY_FCT,
		 const double &TERMIN_GAP_OPENING_PENALTY,
		 const double &TERMIN_GAP_EXTENSION_PENALTY,
		 const Sequence &FIRST_SEQUENCE,
		 const Sequence &SECOND_SEQUENCE,
		 const ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORE,
		 Matrix< DynamicProgrammingMatrixElement> &MATRIX
 )
:
	  m_GapOpeningPenaltyFunction( GAP_OPENING_PENALTY_FCT),   // build member: call constructor of double and initialize with value GAP_OPENING_PENALTY
	  m_GapExtensionPenaltyFunction( GAP_EXTENSION_PENALTY_FCT),
	  m_TerminiGapOpeningPenalty( TERMIN_GAP_OPENING_PENALTY),
	  m_TerminiGapExtensionPenalty( TERMIN_GAP_EXTENSION_PENALTY),
	  m_FirstSequence( FIRST_SEQUENCE),
	  m_SecondSequence( SECOND_SEQUENCE),
	  m_Score( SCORE),
	  m_Matrix( MATRIX)
{}


void NeedlemanWunschAffineGaps::CalculateMatrix()
{
	double
	  minvalue = double ( INT_MIN); // 100000.0 - std::numeric_limits<double>::max();  //  arbitrary choice, negative value mean penalty, first number makes sure that if summed to negative values that it does not get positive (large value)

//	std::cout << __FUNCTION__ << " minvalue: " << minvalue << std::endl;

	DynamicProgrammingMatrixElement
		element,
		* element_ptr,
		element_1,
		element_2,
		element_3;
	double
		max_value,
		local_score;
	std::vector< double>
		profiles;

//	for( size_t i = 0; i < 9; ++i)
//	{
//		m_Matrix( 0, 0).AddValue(i, 0);
//	}

	// each dynamic-programing-matrix-element has a 9 element affine-path-way member
	// first three represent the previous diagonal (times its three possible previous elements)
	// second three go horizontal previous element, last three vertical previous element
	// each triplet follows the same structure


	// HORIZONTAL = columns // nr rows - 1 = FIRST_SEQUENCE.size()
	for( size_t i = 1; i < m_Matrix.GetNumberOfRows(); ++i)
	{
		for( size_t k = 0 ; k <4; ++k)
		{
		  m_Matrix( i, 0).AddValue(k, minvalue);
		}

		//		for( size_t k = 4; k < 5; ++k)
		m_Matrix( i, 0).AddValue(4, -m_TerminiGapOpeningPenalty - double(i - 1) * m_TerminiGapExtensionPenalty);

		for( size_t k = 5; k < 9; ++k)
		m_Matrix( i, 0).AddValue(k, minvalue);
	}

	// VERTICAL = rows // nr cols - 1 = SECOND_SEQUENCE.size()
	for( size_t j = 1; j < m_Matrix.GetNumberOfColumns(); ++j)
	{
		for( size_t k = 0; k < 8; ++k)
		m_Matrix( 0, j).AddValue(k, minvalue);

		//		for( size_t k = 8; k < 9; ++k)
		m_Matrix( 0, j).AddValue(8, -m_TerminiGapOpeningPenalty - double(j - 1) * m_TerminiGapExtensionPenalty);
	}


	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() -1; ++i)  // first sequence
		for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() -1; ++j)  // second sequence
		{
			// DIAGONAL
			element = m_Matrix( i -1, j -1); // access previous diagonal (aligned) element
			local_score =  m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[ j - 1])); // similarity score, note that i is matrix index, matrix has zero in first element, sequence does not, that's why i-1 and j-1 are used here
			m_Matrix( i, j).AddValue( 0, element.BestSubPathWayScore(0) + local_score); // fill array containing affine path ways, first three values from diagonal step
			m_Matrix( i, j).AddValue( 1, element.BestSubPathWayScore(1) + local_score);
			m_Matrix( i, j).AddValue( 2, element.BestSubPathWayScore(2) + local_score);

			// HORIZONTAL
			element = m_Matrix( i -1, j);
			max_value = element.BestSubPathWayScore( 0);
			//profiles = m_FirstSequence[ i - 1].GetProfiles();    // considers only left hand value for threshold, TODO: average with m_SecondSequence[j+1].GetProfiles()
			//profiles = 0.5 * ( m_SecondSequence[ j].GetProfiles() + m_SecondSequence[ j - 1].GetProfiles());
			profiles = MinPerElement( m_SecondSequence[ j].GetProfiles(),m_SecondSequence[ j - 1].GetProfiles());

			m_Matrix( i, j).AddValue( 3, max_value - m_GapOpeningPenaltyFunction->operator ()( profiles));
			max_value = element.BestSubPathWayScore( 1);
			m_Matrix( i, j).AddValue( 4, max_value - m_GapExtensionPenaltyFunction->operator()( profiles));
			max_value = element.BestSubPathWayScore( 2);
			m_Matrix( i, j).AddValue( 5, max_value - m_GapOpeningPenaltyFunction->operator ()( profiles));

			// VERTICAL
			element = m_Matrix( i, j -1);
			max_value = element.BestSubPathWayScore( 0);
			profiles = MinPerElement( m_FirstSequence[ i].GetProfiles(), m_FirstSequence[ i - 1].GetProfiles());
			//profiles = 0.5 * ( m_FirstSequence[ i].GetProfiles() +  m_FirstSequence[ i - 1].GetProfiles());
			//profiles = m_SecondSequence[ j - 1].GetProfiles();
			m_Matrix( i, j).AddValue( 6, max_value -  m_GapOpeningPenaltyFunction->operator ()( profiles));
			max_value = element.BestSubPathWayScore( 1);
			m_Matrix( i, j).AddValue( 7, max_value -  m_GapOpeningPenaltyFunction->operator ()( profiles));
			max_value = element.BestSubPathWayScore( 2);
			m_Matrix( i, j).AddValue( 8, max_value -  m_GapExtensionPenaltyFunction->operator()( profiles));
		}


	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() - 1; ++i)
	{
		element_ptr   = &m_Matrix( i, m_Matrix.GetNumberOfColumns() - 1);
		element_1 = m_Matrix( i - 1, m_Matrix.GetNumberOfColumns() - 1 - 1);
		element_2 = m_Matrix( i - 1, m_Matrix.GetNumberOfColumns() - 1);
		element_3 = m_Matrix( i, m_Matrix.GetNumberOfColumns() - 1 - 1);
		local_score = m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[  m_Matrix.GetNumberOfColumns()-1 - 1]));

		element_ptr->AddValue( 0, element_1.BestSubPathWayScore(0) + local_score);
		element_ptr->AddValue( 1, element_1.BestSubPathWayScore(1) + local_score);
		element_ptr->AddValue( 2, element_1.BestSubPathWayScore(2) + local_score);
		element_ptr->AddValue( 3, element_2.BestSubPathWayScore(0) - m_TerminiGapOpeningPenalty);
		element_ptr->AddValue( 4, element_2.BestSubPathWayScore(1) - m_TerminiGapExtensionPenalty);
		element_ptr->AddValue( 5, element_2.BestSubPathWayScore(2) - m_TerminiGapOpeningPenalty);

		profiles = MinPerElement(m_FirstSequence[i].GetProfiles(),m_FirstSequence[i-1].GetProfiles() );

		//profiles = 0.5 * ( m_FirstSequence[ i].GetProfiles() +  m_FirstSequence[ i - 1].GetProfiles()); // TODO : CHECK
		element_ptr->AddValue( 6, element_3.BestSubPathWayScore(0) - ( *m_GapOpeningPenaltyFunction)( profiles));
		element_ptr->AddValue( 7, element_3.BestSubPathWayScore(1) - ( *m_GapOpeningPenaltyFunction)( profiles));
		element_ptr->AddValue( 8, element_3.BestSubPathWayScore(2) - ( *m_GapExtensionPenaltyFunction)( profiles));
	}


	for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() - 1; ++j)
	{
		element_ptr   = &m_Matrix( m_Matrix.GetNumberOfRows() - 1, j);
		element_1 = m_Matrix( m_Matrix.GetNumberOfRows() - 1 - 1, j - 1);
		element_2 = m_Matrix( m_Matrix.GetNumberOfRows() - 1 - 1, j);
		element_3 = m_Matrix( m_Matrix.GetNumberOfRows() - 1, j - 1);
		local_score = m_Score->operator()( std::make_pair( m_FirstSequence[  m_Matrix.GetNumberOfRows()-1 - 1], m_SecondSequence[ j - 1]));

		element_ptr->AddValue( 0, element_1.BestSubPathWayScore(0) + local_score);
		element_ptr->AddValue( 1, element_1.BestSubPathWayScore(1) + local_score);
		element_ptr->AddValue( 2, element_1.BestSubPathWayScore(2) + local_score);

		profiles = MinPerElement(m_SecondSequence[ j].GetProfiles(),m_SecondSequence[ j - 1].GetProfiles());
		//profiles = 0.5 * ( m_SecondSequence[ j].GetProfiles() + m_SecondSequence[ j - 1].GetProfiles()); // TODO : CHECK
		element_ptr->AddValue( 3, element_2.BestSubPathWayScore(0) - m_GapOpeningPenaltyFunction->operator ()( profiles));
		element_ptr->AddValue( 4, element_2.BestSubPathWayScore(1) - m_GapExtensionPenaltyFunction->operator ()( profiles));
		element_ptr->AddValue( 5, element_2.BestSubPathWayScore(2) - m_GapOpeningPenaltyFunction->operator ()( profiles));
		element_ptr->AddValue( 6, element_3.BestSubPathWayScore(0) - m_TerminiGapOpeningPenalty);
		element_ptr->AddValue( 7, element_3.BestSubPathWayScore(1) - m_TerminiGapOpeningPenalty);
		element_ptr->AddValue( 8, element_3.BestSubPathWayScore(2) - m_TerminiGapExtensionPenalty);
	}

	element_ptr   = &m_Matrix( m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns() -1);
	element_1     = m_Matrix( m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1-1);
	element_2 	  = m_Matrix( m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1);
	element_3     = m_Matrix( m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns()-1-1);
	local_score   = m_Score->operator()( std::make_pair( m_FirstSequence[  m_Matrix.GetNumberOfRows()-1 - 1], m_SecondSequence[ m_Matrix.GetNumberOfColumns()-1 - 1]));

	element_ptr->AddValue(0, element_1.BestSubPathWayScore(0)+ local_score);
	element_ptr->AddValue(1, element_1.BestSubPathWayScore(1)+ local_score);
	element_ptr->AddValue(2, element_1.BestSubPathWayScore(2)+ local_score);
	element_ptr->AddValue(3, element_2.BestSubPathWayScore(0)- m_TerminiGapOpeningPenalty);
	element_ptr->AddValue(4, element_2.BestSubPathWayScore(1)- m_TerminiGapExtensionPenalty);
	element_ptr->AddValue(5, element_2.BestSubPathWayScore(2)- m_TerminiGapOpeningPenalty);
	element_ptr->AddValue(6, element_3.BestSubPathWayScore(0)- m_TerminiGapOpeningPenalty);
	element_ptr->AddValue(7, element_3.BestSubPathWayScore(1)- m_TerminiGapOpeningPenalty);
	element_ptr->AddValue(8, element_3.BestSubPathWayScore(2)- m_TerminiGapExtensionPenalty);
}

std::pair< double, std::vector< std::pair< int, int> > >
NeedlemanWunschAffineGaps:: TraceBack() const
{

	double
		best = -std::numeric_limits< double>::max(),
		best2;

	int
		max = std::numeric_limits<int>::max(),
		position = 0,
		group,
		start,
		end;

	std::list<std::pair<int, int> >
		alignment;


	{
		DynamicProgrammingMatrixElement
			element = m_Matrix(m_Matrix.GetNumberOfRows() - 1, m_Matrix.GetNumberOfColumns() - 1);

		// get best scoring last matrix element
		for (size_t i = 0; i < 9; i++)
		{
			if ( element.GetAffinePathWay( i) > best)
			{
				best = element.GetAffinePathWay( i);
				position = i;
			}
		}
	}


	size_t i = m_Matrix.GetNumberOfRows()    - 1;
	size_t j = m_Matrix.GetNumberOfColumns() - 1;
	group =  int (double (position) / 3.0);
	position -= 3 * group;

//	std::cout << __FUNCTION__ << ": last element: " << m_Matrix( m_Matrix.GetNumberOfRows() - 1, m_Matrix.GetNumberOfColumns() - 1) << "\n";
//
//	std::cout << "group: " << group << " position: " << position << std::endl;

	while (i > 0 || j > 0)
	{
//	  std::cout << __FUNCTION__ << ": group: " << group << " (position: " << position << ")" << std::endl;
		if (group == 0)
		{
			alignment.push_front(std::make_pair(i - 1, j - 1));
			--i;
			--j;
		}

		else if (group == 1)
		{
			alignment.push_front(std::make_pair(i - 1, max));
			--i;
		}

		else if (group == 2)
		{
			alignment.push_front(std::make_pair(max, j - 1));
			--j;
		}

		group = position;
		start = position * 3;
		end = (position * 3) + 3;
//		std::cout << "group: " << group << " position: " << position << " start: " << start << " end: " << end << " i: " << i << " j: " << j << std::endl;

		best2 = -std::numeric_limits< double>::max();

		for (int count = start; count < end; ++count)
		{
			if (m_Matrix(i, j).GetAffinePathWay( count) > best2)
			{
				best2 = m_Matrix(i, j).GetAffinePathWay( count);
				position = count;
			}
		}
//		std::cout << "pos: " << position << std::endl;
		position -= start;
//		std::cout << "position: " << position << std::endl << std::endl;
	}

//	std::cout << __FUNCTION__ << ": last element: " << m_Matrix( m_Matrix.GetNumberOfRows() - 1, m_Matrix.GetNumberOfColumns() - 1) << "\n";

	std::vector< std::pair< int, int> >
		converted( alignment.begin(), alignment.end());
	return std::make_pair( best, converted);  // return last matrix element as first of pair
}


const Matrix< DynamicProgrammingMatrixElement> &
NeedlemanWunschAffineGaps::GetMatrix() const
{
	return m_Matrix;
}


