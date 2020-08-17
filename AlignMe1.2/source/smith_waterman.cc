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
//!  Smith Waterman algorithm.
//! As described in "Biological sequence analysis"
//! from Durbin, Eddy, Krogh, Mitchison
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#include "../include/smith_waterman.h"


// default constructor
SmithWaterman::SmithWaterman()
  : m_Matrix( 0, 0)
{}


// construct from data
SmithWaterman::SmithWaterman
( 
		const double &GAP_OPENING_PENALTY,
		const double &GAP_EXTENSION_PENALTY,
		const Sequence &FIRST_SEQUENCE,
		const Sequence &SECOND_SEQUENCE,
		const ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORE,
		Matrix< DynamicProgrammingMatrixElement> &MATRIX
)
	:
	m_GapOpeningPenalty( GAP_OPENING_PENALTY),   // build member: call constructor of double and initialize with value GAP_OPENING_PENALTY
	m_GapExtensionPenalty( GAP_EXTENSION_PENALTY),
	m_FirstSequence( FIRST_SEQUENCE),
	m_SecondSequence( SECOND_SEQUENCE),
	m_Score( SCORE),
	m_Matrix( MATRIX),
	m_BestScoringElement(0, 0)
	{}

// copy constructor


void SmithWaterman::CalculateMatrix()
{
	// first column
	for( size_t i = 0; i < m_Matrix.GetNumberOfRows(); ++i)
	{
		m_Matrix( i, 0) = DynamicProgrammingMatrixElement( 0.0, std::numeric_limits< size_t>::max(), std::numeric_limits< size_t>::max());
		DebugWrite( i <<  " 0  " <<  m_Matrix( i, 0));
	}

	// fill first line, except first element
	for( size_t j = 1; j < m_Matrix.GetNumberOfColumns(); ++j)
	{
		m_Matrix( 0, j) = DynamicProgrammingMatrixElement( 0.0, std::numeric_limits< size_t>::max(), std::numeric_limits< size_t>::max());
		DebugWrite( "0  " << j << " 0  " <<  m_Matrix( 0, j));
	}
  
	// calculate matrix without
	for( size_t i = 1; i < m_Matrix.GetNumberOfRows(); ++i)
		for( size_t j = 1; j < m_Matrix.GetNumberOfColumns(); ++j)
		{
			m_Matrix( i, j) = FindBestScoringPriorNeighbor( i, j);
			DebugWrite( i << " " << j << "   " << m_Matrix( i, j));
		}

}


DynamicProgrammingMatrixElement
SmithWaterman::FindBestScoringPriorNeighbor( const size_t &I, const size_t &J)
{
	double
		penalty;

	DynamicProgrammingMatrixElement
		element,
		best_element;

	// evaluate diagonal back element
	best_element = DynamicProgrammingMatrixElement(  m_Matrix( I - 1, J - 1).GetValue() + m_Score->operator()( std::make_pair( m_FirstSequence[ I - 1], m_SecondSequence[ J - 1])), I - 1, J - 1);

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
  
	//DebugWrite ( "element: " << element);
	// check for best scoring element
	if( element.GetValue() > best_element.GetValue())
	{
		best_element = element;
	}
	//DebugWrite ( "best element: " <<  best_element);
	if( best_element.GetValue() <= 0) // equal zero: no pointer to
	{
		best_element = DynamicProgrammingMatrixElement( 0.0, std::numeric_limits< size_t>::max(), std::numeric_limits< size_t>::max());
		//DebugWrite ( "best element (was negative) is set to: " <<  best_element);
	}
	else if( m_Matrix( m_BestScoringElement.first, m_BestScoringElement.second).GetValue() < best_element.GetValue())
	{
		//DebugWrite ("m_Matrix: " << m_Matrix( m_BestScoringElement.first, m_BestScoringElement.second).GetValue());
		//DebugWrite ( "best_element.GetValue(): " << best_element.GetValue());
		m_BestScoringElement = std::make_pair( I, J);
		DebugWrite ( "best element (was positive) is set to position: " << I << " " << J);
	}

	// assign element
	return best_element;
}

			 

std::pair< double, std::vector< std::pair< int, int> > >
SmithWaterman::TraceBack() const
{
	size_t
		i( m_BestScoringElement.first),
		j( m_BestScoringElement.second);

	int
		max( std::numeric_limits< int>::max());
  
	std::list< std::pair< int, int> >
		alignment;
	//DebugWrite ( "best index i: " << i);
	//DebugWrite ( "best index j: " << j);
	if( i == std::numeric_limits< size_t>::max() || j ==  std::numeric_limits< size_t>::max())
	{
		return std::make_pair( 0.0, std::vector< std::pair< int, int> >());
	}
  
	while( m_Matrix( i, j).GetValue() > 0.0)
	{
		if( m_Matrix( i, j).GetIndicesOfPreviousElement().first < i && m_Matrix( i, j).GetIndicesOfPreviousElement().second < j)
		{
			DebugWrite( i - 1 << " " << j - 1);
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
SmithWaterman::GetMatrix() const
{
	return m_Matrix;
}


