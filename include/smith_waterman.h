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
//!  The algorithm for local alignments.
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

#include <list>
#include <limits>

#include <boost/shared_ptr.hpp>

#include "sequence.h"
#include "dynamic_programing_matrix_element.h"
#include "function.t.h"
#include "matrix.t.h"
#include "macro_functions_read_write.h"

class SmithWaterman
{
private:
	double                                            m_GapOpeningPenalty;
	double                                            m_GapExtensionPenalty;
	Sequence                                        m_FirstSequence;
	Sequence                                        m_SecondSequence;
	ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >   m_Score;
	Matrix< DynamicProgrammingMatrixElement>          m_Matrix;
	std::pair< size_t, size_t>                        m_BestScoringElement;

public:

	// default constructor
	SmithWaterman();

	// construct from data
	SmithWaterman
	(
			const double &GAP_OPENING_PENALTY,
			const double &GAP_EXTENSION_PENALTY,
			const Sequence &FIRST_SEQUENCE,
			const Sequence &SECOND_SEQUENCE,
			const ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORE,
			Matrix< DynamicProgrammingMatrixElement> &MATRIX
	);

    // destructor
	~SmithWaterman(){}

	//! calculates all possible pathways (alignments)
	void CalculateMatrix();
    
	//! identifies the actual alignment (highest scoring pathway) after matrix has been calculated
	std::pair< double, std::vector< std::pair< int, int> > >
	TraceBack() const;

private:
	DynamicProgrammingMatrixElement
	FindBestScoringPriorNeighbor( const size_t &I, const size_t &J);

public:
	Matrix< DynamicProgrammingMatrixElement>
	GetMatrix() const;


    std::ostream &
    Write( std::ostream &STREAM) const
    {
    	//	int value;
    	STREAM << "Smith-Waterman-Local-Alignment" << "\n";
    	STREAM << "gap_opening_penalty: " <<  m_GapOpeningPenalty << "\n";
    	STREAM << "gap_extension_penalty: " << m_GapExtensionPenalty << "\n";
    	STREAM << "best_scoring_element: " << m_BestScoringElement.first << " " << m_BestScoringElement.second << "\n";
    	for( size_t i = 0; i < m_Matrix.GetNumberOfRows(); ++i)
    	{
    		for( size_t j = 0; j < m_Matrix.GetNumberOfColumns(); ++j)
    		{
    			DynamicProgrammingMatrixElement element( m_Matrix( i, j));
    			STREAM.width( 5);
    			STREAM << element.GetValue() << " (";
    			if( element.GetIndicesOfPreviousElement().first == std::numeric_limits< size_t>::max())
    			{
    				STREAM.width( 3);
    				STREAM << "x";
    			}
    			else
    			{
    				STREAM.width( 3);
    				STREAM << element.GetIndicesOfPreviousElement().first;
    			}
    			STREAM << ",";
    			if( element.GetIndicesOfPreviousElement().second == std::numeric_limits< size_t>::max())
    			{
    				STREAM.width( 3);
    				STREAM << "x";
    			}
    			else
    			{
    				STREAM.width( 3);
    				STREAM << element.GetIndicesOfPreviousElement().second;
    			}
    			STREAM <<  ")     ";
    		}
    		STREAM << "\n";
    	}
    	return STREAM;
    }

};

#endif
