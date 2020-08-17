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
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef NEEDLMAN_WUNSCH_AFFINE_GAPS_H
#define NEEDLMAN_WUNSCH_AFFINE_GAPS_H

#include "sequence.h"
#include "dynamic_programing_matrix_element.h"
#include "function.t.h"
#include "matrix.t.h"
#include "macro_functions_read_write.h"
#include <boost/shared_ptr.hpp>

#include <list>

class NeedlemanWunschAffineGaps
{
private:
	ShPtr< Function< std::vector< double>, double> >   m_GapOpeningPenaltyFunction;
	ShPtr< Function< std::vector< double>, double> >   m_GapExtensionPenaltyFunction;
	double                                                         m_TerminiGapOpeningPenalty;
	double                                                         m_TerminiGapExtensionPenalty;
	Sequence                                                     m_FirstSequence;
	Sequence                                                     m_SecondSequence;
	ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >   m_Score;
	Matrix< DynamicProgrammingMatrixElement>                       m_Matrix;

public:

	// default constructor
	NeedlemanWunschAffineGaps();

	// construct from data
	NeedlemanWunschAffineGaps
	(
			const ShPtr< Function< std::vector< double>, double> > &GAP_OPENING_PENALTY_FCT,
			const ShPtr< Function< std::vector< double>, double> > &GAP_EXTENSION_PENALTY_FCT,
			const double &TERMIN_GAP_OPENING_PENALTY,
			const double &TERMIN_GAP_EXTENSION_PENALTY,
			const Sequence &FIRST_SEQUENCE,
			const Sequence &SECOND_SEQUENCE,
			const ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORE,
			Matrix< DynamicProgrammingMatrixElement> &MATRIX
	);

	// copy constructor

    // destructor
	~NeedlemanWunschAffineGaps(){}

    
	void CalculateMatrix();
    
	std::pair< double, std::vector< std::pair< int, int> > >
	TraceBack() const;

private:
	DynamicProgrammingMatrixElement
	FindBestScoringPriorNeighbor( const size_t &I, const size_t &J);

	DynamicProgrammingMatrixElement
	FindBestScoringPriorNeighborInLastRow( const size_t &J);

	DynamicProgrammingMatrixElement
	FindBestScoringPriorNeighborInLastColumn( const size_t &I);

//	DynamicProgrammingMatrixElement
//	FindBestScoringPriorNeighborWithTerminiGapPenalties( const size_t &I, const size_t &J);

public:
	const Matrix< DynamicProgrammingMatrixElement> &
	GetMatrix() const;


//    std::ostream &
//    Write( std::ostream &STREAM) const
//    {
//    	//	int value;
//    	STREAM << "Needleman-Wunsch-Global-Alignment" << "\n";
//    	STREAM << "gap_opening_penalty: " <<  m_GapOpeningPenalty << "\n";
//    	STREAM << "gap_extension_penalty: " << m_GapExtensionPenalty << "\n";
//    	for( size_t i = 0; i < m_Matrix.GetNumberOfRows(); ++i)
//    	{
//    		for( size_t j = 0; j < m_Matrix.GetNumberOfColumns(); ++j)
//    		{
//    			DynamicProgrammingMatrixElement element( m_Matrix( i, j));
//    			STREAM.width( 5);
//    			STREAM << element.GetValue() << " (";
//    			if( element.GetIndicesOfPreviousElement().first == std::numeric_limits< size_t>::max())
//    			{
//    				STREAM.width( 3);
//    				STREAM << "x";
//    			}
//    			else
//    			{
//    				STREAM.width( 3);
//    				STREAM << element.GetIndicesOfPreviousElement().first;
//    			}
//    			STREAM << ",";
//    			if( element.GetIndicesOfPreviousElement().second == std::numeric_limits< size_t>::max())
//    			{
//    				STREAM.width( 3);
//    				STREAM << "x";
//    			}
//    			else
//    			{
//    				STREAM.width( 3);
//    				STREAM << element.GetIndicesOfPreviousElement().second;
//    			}
//    			STREAM <<  ")     ";
//    		}
//    		STREAM << "\n";
//    	}
//    	return STREAM;
//    }

};

#endif
