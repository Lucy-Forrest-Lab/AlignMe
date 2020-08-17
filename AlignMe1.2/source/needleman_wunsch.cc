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
//!  Needleman Wunsch algorithm.
//! As described in "Biological sequence analysis"
//! from Durbin, Eddy, Krogh, Mitchison
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#include "../include/needleman_wunsch.h"
#include <limits>

// default constructor
NeedlemanWunsch::NeedlemanWunsch()
  : m_Matrix( 0, 0)
{}


// construct from data
// needleman_wunsch( gap_opening_penalty, gap_extension_penalty, termini_gap_opening_penalty, termini_gap_extension_penalty,
// first_sequence, second_sequence, scores);
// assign these values directly to classmembers with an initialization-list

NeedlemanWunsch::NeedlemanWunsch
( 
 const double &GAP_OPENING_PENALTY, 
 const double &GAP_EXTENSION_PENALTY,
 const double &TERMIN_GAP_OPENING_PENALTY,
 const double &TERMIN_GAP_EXTENSION_PENALTY,
 const Sequence &FIRST_SEQUENCE,
 const Sequence &SECOND_SEQUENCE,
 const ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORE,
 Matrix< DynamicProgrammingMatrixElement> &MATRIX
 )
  :
  m_GapOpeningPenalty( GAP_OPENING_PENALTY),   // build member: call constructor of double and initialize with value GAP_OPENING_PENALTY
  m_GapExtensionPenalty( GAP_EXTENSION_PENALTY),
  m_TerminiGapOpeningPenalty( TERMIN_GAP_OPENING_PENALTY),
  m_TerminiGapExtensionPenalty( TERMIN_GAP_EXTENSION_PENALTY),
  m_FirstSequence( FIRST_SEQUENCE),
  m_SecondSequence( SECOND_SEQUENCE),
  m_Score( SCORE),
  m_Matrix( MATRIX)
{}


// copy constructor


void NeedlemanWunsch::CalculateMatrix()
{
  
	// first element
	m_Matrix( 0, 0) = DynamicProgrammingMatrixElement( 0.0, std::numeric_limits< size_t>::max(), std::numeric_limits< size_t>::max());
  
	// gap opening penalty for beginning for first sequence
	m_Matrix( 1, 0) = DynamicProgrammingMatrixElement( -1.0 * m_TerminiGapOpeningPenalty, 0, 0);

	// gap opening penalty for beginning for second sequence
	m_Matrix( 0, 1) = DynamicProgrammingMatrixElement( -1.0 * m_TerminiGapOpeningPenalty, 0, 0);

	// first column, except first element
	for( size_t i = 2; i < m_Matrix.GetNumberOfRows(); ++i)
	{
		m_Matrix( i, 0) = DynamicProgrammingMatrixElement( m_Matrix( 1, 0).GetValue() - 1.0 * (i - 1) * m_TerminiGapExtensionPenalty , i - 1, 0);
		DebugWrite( i << " 0  " << -1.0 * i * m_TerminiGapExtensionPenalty << "  "<< m_Matrix( i, 0));
	}

	// fill first line, except first element
	for( size_t j = 2; j < m_Matrix.GetNumberOfColumns(); ++j)
	{
		m_Matrix( 0, j) = DynamicProgrammingMatrixElement( m_Matrix( 0, 1).GetValue() - 1.0 * (j - 1) * m_TerminiGapExtensionPenalty, 0 , j - 1);
		DebugWrite( "0 " << j << "   "  << -1.0 * j * m_TerminiGapExtensionPenalty << "  "<<  m_Matrix( 0, j));
	}
  

	// calculate matrix without last
	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() - 1; ++i)
		for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() - 1; ++j)
		{
			m_Matrix( i, j) = FindBestScoringPriorNeighbor( i, j);
			DebugWrite( i << " " << j << "   " << m_Matrix( i, j));
		}

	size_t
		last_row( m_Matrix.GetNumberOfRows() - 1),
		last_col( m_Matrix.GetNumberOfColumns() - 1);


		// calculate matrix
		for( size_t i = 1; i < m_Matrix.GetNumberOfRows() - 1; ++i)
		{
//			m_Matrix( i, last_col) = FindBestScoringPriorNeighborWithTerminiGapPenalties( i, last_col);
			m_Matrix( i, last_col) = FindBestScoringPriorNeighborInLastColumn( i);
			DebugWrite( " last col wtgp: " << i << " " << last_col << "   " << m_Matrix( i, last_col));
		}

		for( size_t j = 1; j < m_Matrix.GetNumberOfColumns(); ++j)
		{
			m_Matrix( last_row, j) = FindBestScoringPriorNeighborInLastRow( j);
			DebugWrite( "last row wtgp: " << last_row << " " << j << "   " << m_Matrix( last_row, j));
		}
//	}
}

DynamicProgrammingMatrixElement
NeedlemanWunsch::FindBestScoringPriorNeighbor( const size_t &I, const size_t &J)
{
	double
		penalty;

	DynamicProgrammingMatrixElement
		element,
		best_element;
  
	// evaluate diagonal back element, sum its value with the similarity score at this position
	best_element = DynamicProgrammingMatrixElement(  m_Matrix( I - 1, J - 1).GetValue() + m_Score->operator()( std::make_pair( m_FirstSequence[ I - 1], m_SecondSequence[ J - 1])), I - 1, J - 1); // sequence does start at index i = 1, j = 1 respectively, therefore Sequence[ I - 1]
	DebugWrite( "diagonal value:" << best_element.GetValue() << " = " << m_Matrix( I - 1, J - 1).GetValue() << " + " << m_Score->operator()( std::make_pair( m_FirstSequence[ I - 1], m_SecondSequence[ J - 1])));
	DebugWrite( "prev: " << best_element.GetIndicesOfPreviousElement());
	// evaluate upwards back element
	// test whether gap opening or extension penalty
	if( m_Matrix( I - 1, J).GetIndicesOfPreviousElement().second == J)
	{
		penalty = m_GapExtensionPenalty;
	}
	else
	{
		penalty = m_GapOpeningPenalty;
	}

	element = DynamicProgrammingMatrixElement(  m_Matrix( I - 1, J).GetValue() - penalty, I - 1, J);
	DebugWrite( "upwards value: " << element.GetValue() << " = " << m_Matrix( I - 1, J).GetValue() << " - " << penalty);
	DebugWrite( "prev: " << element.GetIndicesOfPreviousElement());

	// check for best scoring element
	if( element.GetValue() > best_element.GetValue())
	{
		best_element = element;
	}

	// evaluate left back element
	if( m_Matrix( I, J - 1).GetIndicesOfPreviousElement().first == I)
	{
		penalty = m_GapExtensionPenalty;
	}
	else
	{
		penalty = m_GapOpeningPenalty;
	}

	element = DynamicProgrammingMatrixElement(  m_Matrix( I, J - 1).GetValue() - penalty, I, J - 1);

	DebugWrite( "left value: " << element.GetValue() << " = " << m_Matrix( I, J - 1).GetValue() << " - " << penalty);
	// check for best scoring element
	if( element.GetValue() > best_element.GetValue())
	{
		best_element = element;
	}

	DebugWrite( "best element " << best_element);

	// assign element
	return best_element;
}

DynamicProgrammingMatrixElement
NeedlemanWunsch::FindBestScoringPriorNeighborInLastColumn( const size_t &I)
{
	double
		penalty;

	DynamicProgrammingMatrixElement
		element,
		best_element;

	size_t
		j( m_Matrix.GetNumberOfColumns() - 1);

	// evaluate diagonal back element, sum its value with the similarity score at this position
	best_element = DynamicProgrammingMatrixElement(  m_Matrix( I - 1, j - 1).GetValue() + m_Score->operator()( std::make_pair( m_FirstSequence[ I - 1], m_SecondSequence[ j - 1])), I - 1, j - 1); // sequence does start at index i = 1, j = 1 respectively, therefore Sequence[ I - 1]

	DebugWrite( "diagonal value:" << best_element.GetValue() << " = " << m_Matrix( I - 1, j - 1).GetValue() << " + " << m_Score->operator()( std::make_pair( m_FirstSequence[ I - 1], m_SecondSequence[ j - 1])));


	// evaluate upwards back element
	// test whether gap opening or extension penalty
	if( m_Matrix( I - 1, j).GetIndicesOfPreviousElement().second == j)
	{
		penalty = m_TerminiGapExtensionPenalty;
	}
	else
	{
		penalty = m_TerminiGapOpeningPenalty;
	}

	element = DynamicProgrammingMatrixElement(  m_Matrix( I - 1, j).GetValue() - penalty, I - 1, j);

	DebugWrite( "upwards value: " << element.GetValue() << " = " << m_Matrix( I - 1, j).GetValue() << " - " << penalty);

	// check for best scoring element
	if( element.GetValue() > best_element.GetValue())
	{
		DebugWrite( "upward gap better than diagonal");
		best_element = element;
	}

	// evaluate left back element
	if( m_Matrix( I, j - 1).GetIndicesOfPreviousElement().first == I)
	{
		penalty = m_GapExtensionPenalty;
		DebugWrite( "why a gap extension here in the last col?");
	}
	else
	{
		penalty = m_GapOpeningPenalty;
		DebugWrite( "why a gap opening here in the last col?");
	}

	element = DynamicProgrammingMatrixElement(  m_Matrix( I, j - 1).GetValue() - penalty, I, j - 1);

	DebugWrite( "left value: " << element.GetValue() << " = " << m_Matrix( I, j - 1).GetValue() << " - " << penalty);

	// check for best scoring element
	if( element.GetValue() > best_element.GetValue())
	{
		DebugWrite( "why did he win in the last col???");
		best_element = element;
	}

	// assign element
	return best_element;
}



DynamicProgrammingMatrixElement
NeedlemanWunsch::FindBestScoringPriorNeighborInLastRow( const size_t &J)
{
	double
		penalty;

	DynamicProgrammingMatrixElement
		element,
		best_element;

	size_t i( m_Matrix.GetNumberOfRows() - 1);

	// evaluate diagonal back element, sum its value with the similarity score at this position
	best_element = DynamicProgrammingMatrixElement(  m_Matrix( i - 1, J - 1).GetValue() + m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[ J - 1])), i - 1, J - 1); // sequence does start at index i = 1, j = 1 respectively, therefore Sequence[ i - 1]

	DebugWrite( "diagonal value:" << best_element.GetValue() << " = " << m_Matrix( i - 1, J - 1).GetValue() << " + " << m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[ J - 1])));


	// evaluate upwards back element
	// test whether gap opening or extension penalty
	if( m_Matrix( i - 1, J).GetIndicesOfPreviousElement().second == J)
	{
		DebugWrite( "why a gap extension here in the last row?");
		penalty = m_GapExtensionPenalty;
	}
	else
	{
		DebugWrite( "why a gap opening here in the last row?");
		penalty = m_GapOpeningPenalty;
	}

	element = DynamicProgrammingMatrixElement(  m_Matrix( i - 1, J).GetValue() - penalty, i - 1, J);

	DebugWrite( "upwards value: " << element.GetValue() << " = " << m_Matrix( i - 1, J).GetValue() << " - " << penalty);

	// check for best scoring element
	if( element.GetValue() > best_element.GetValue())
	{
		DebugWrite( "why upward gap better than diagonal???");
		best_element = element;
	}

	// evaluate left back element
	if( m_Matrix( i, J - 1).GetIndicesOfPreviousElement().first == i)
	{
		penalty = m_TerminiGapExtensionPenalty;
	}
	else
	{
		penalty = m_TerminiGapOpeningPenalty;
	}

	element = DynamicProgrammingMatrixElement(  m_Matrix( i, J - 1).GetValue() - penalty, i, J - 1);

	DebugWrite( "left value: " << element.GetValue() << " = " << m_Matrix( i, J - 1).GetValue() << " - " << penalty);

	// check for best scoring element
	if( element.GetValue() > best_element.GetValue())
	{
		DebugWrite( "left gap value wins");
		best_element = element;
	}

	// assign element
	return best_element;
}



std::pair< double, std::vector< std::pair< int, int> > >
NeedlemanWunsch:: TraceBack() const
{
	size_t
		i( m_Matrix.GetNumberOfRows() - 1),
		j( m_Matrix.GetNumberOfColumns() - 1);

	int
		max( std::numeric_limits< int>::max());
  
	std::list< std::pair< int, int> >
		alignment;
  
  
	while( i > 0 || j > 0) // omit
	{
		if( m_Matrix( i, j).GetIndicesOfPreviousElement().first < i && m_Matrix( i, j).GetIndicesOfPreviousElement().second < j)
		{
			DebugWrite( "diagonal " << i - 1 << " " << j - 1);
			alignment.push_front( std::make_pair( i - 1, j - 1)); // alignment contains indices of sequence, not of matrix!
			--i;
			--j;
		}
		else if( m_Matrix( i, j).GetIndicesOfPreviousElement().first == i && m_Matrix( i, j).GetIndicesOfPreviousElement().second < j)
		{
			DebugWrite( max << " " << j - 1);
			alignment.push_front( std::make_pair( max, j - 1));   // use max value for gap identification
			--j;
		}
		else if( m_Matrix( i, j).GetIndicesOfPreviousElement().first < i && m_Matrix( i, j).GetIndicesOfPreviousElement().second == j)
		{
			DebugWrite( i - 1 << " " << max);
			alignment.push_front( std::make_pair( i - 1, max));  // use max value for gap identification
			--i;
		}
	}
	std::vector< std::pair< int, int> >
		converted( alignment.begin(), alignment.end());
	return std::make_pair( m_Matrix( m_Matrix.GetNumberOfRows() - 1, m_Matrix.GetNumberOfColumns() - 1).GetValue(), converted);
}


Matrix< DynamicProgrammingMatrixElement>
NeedlemanWunsch::GetMatrix() const
{
	return m_Matrix;
}


