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
		 const boost::shared_ptr< Function< std::vector< double>, double> > &GAP_OPENING_PENALTY_FCT,
		 const boost::shared_ptr< Function< std::vector< double>, double> > &GAP_EXTENSION_PENALTY_FCT,
		 const double &TERMIN_GAP_OPENING_PENALTY,
		 const double &TERMIN_GAP_EXTENSION_PENALTY,
		 const AASequence &FIRST_SEQUENCE,
		 const AASequence &SECOND_SEQUENCE,
		 const boost::shared_ptr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORE,
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
	int
		minvalue = INT_MIN;
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

	for( size_t i = 0; i < 9; ++i)
	{
		m_Matrix( 0, 0).AddValue(i, 0);
	}

	for( size_t i = 1; i < m_Matrix.GetNumberOfRows(); ++i)
	{
		for( size_t j = 0 ; j <4; ++j)
		m_Matrix( i, 0).AddValue(j, minvalue);

		for( size_t j = 4; j < 5; ++j)
		m_Matrix( i, 0).AddValue(j, -1.0 * m_TerminiGapOpeningPenalty - 1.0 * (i - 1) * m_TerminiGapExtensionPenalty);

		for( size_t j = 5; j < 9; ++j)
		m_Matrix( i, 0).AddValue(j, minvalue);
	}

	for( size_t i = 1; i < m_Matrix.GetNumberOfColumns(); ++i)
	{
		for( size_t j = 0; j < 8; ++j)
		m_Matrix( 0, i).AddValue(j, minvalue);

		for( size_t j = 8; j < 9; ++j)
		m_Matrix( 0, i).AddValue(j, -1.0 * m_TerminiGapOpeningPenalty - 1.0 * (i - 1) * m_TerminiGapExtensionPenalty);
	}


	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() -1; ++i)
		for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() -1; ++j)
		{
			element = m_Matrix( i -1, j -1);
			local_score =  m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[ j - 1]));
			m_Matrix( i, j).AddValue( 0, element.BestSubPathWayScore(0) + local_score);
			m_Matrix( i, j).AddValue( 1, element.BestSubPathWayScore(1) + local_score);
			m_Matrix( i, j).AddValue( 2, element.BestSubPathWayScore(2) + local_score);

			element = m_Matrix( i -1, j);
			max_value = element.BestSubPathWayScore( 0);
			//profiles = m_FirstSequence[ i].GetProfiles();
			profiles = m_SecondSequence[ j].GetProfiles();
			m_Matrix( i, j).AddValue( 3, max_value - m_GapOpeningPenaltyFunction->operator ()( profiles));
			max_value = element.BestSubPathWayScore( 1);
			m_Matrix( i, j).AddValue( 4, max_value - m_GapExtensionPenaltyFunction->operator()( profiles));
			max_value = element.BestSubPathWayScore( 2);
			m_Matrix( i, j).AddValue( 5, max_value - m_GapOpeningPenaltyFunction->operator ()( profiles));

			element = m_Matrix( i, j -1);
			max_value = element.BestSubPathWayScore( 0);
			profiles = m_FirstSequence[ i].GetProfiles();
			//profiles = m_SecondSequence[ j].GetProfiles();
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

std::pair< double, std::list< std::pair< int, int> > >
NeedlemanWunschAffineGaps:: TraceBack() const
{

	double
		best = INT_MIN;

	int
		max(std::numeric_limits<int>::max()),
		position,
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
	position =  position - (3 * group);

	while (i > 0 || j > 0)
	{
		if (group == 0)
		{
			alignment.push_front(std::make_pair(i - 1, j - 1));
			i--;
			j--;
		}

		else if (group == 1)
		{
			alignment.push_front(std::make_pair(i - 1, max));
			i--;
		}

		else if (group == 2)
		{
			alignment.push_front(std::make_pair(max, j - 1));
			j--;
		}

		group = position;
		start = position * 3;
		end = (position * 3) + 3;

		double best2 = INT_MIN;

		for (int count = start; count < end; count++)
		{
			if (m_Matrix(i, j).GetAffinePathWay( count) > best2)
			{
				best2 = m_Matrix(i, j).GetAffinePathWay( count);
				position = count;
			}
		}
		position = position - start;
	}
	//std::count << "last element: " << m_Matrix( m_Matrix.GetNumberOfRows() - 1, m_Matrix.GetNumberOfColumns() - 1) << "\n";

	return std::make_pair( best, alignment);  // return last matrix element as first of pair
}


const Matrix< DynamicProgrammingMatrixElement> &
NeedlemanWunschAffineGaps::GetMatrix() const
{
	return m_Matrix;
}


