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
//!  anchor.h
//!
//! Anchors allow to define residues that are kept aligned.
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////



#ifndef ANCHOR_H_
#define ANCHOR_H_

#include "matrix.t.h"
#include "dynamic_programing_matrix_element.h"
#include "triplet.t.h"

std::vector< Triplet< int, int, double> >
BuildAnchorMatrix( std::istream &STREAM, Matrix< DynamicProgrammingMatrixElement> &MATRIX, const double &PENALTY = -1.0)
{
	// anchors need to be given in usual vector notation
	double
		shift_penalty;
	size_t
		nr, ii, jj;
	std::vector< Triplet<int,int,double> >
		anchors;
	STREAM >> nr;
	for( size_t r = 0; r < nr; ++r)
	{
		STREAM >> ii >> jj >> shift_penalty;
		anchors.push_back( Triplet<int,int,double>( ii, jj, shift_penalty));
//		ii -= 1;
//		jj -= 1;
		ii += 1;  // first residue id is 0
		jj += 1;
		shift_penalty *= -1;  // using positive values in file!

		if( ii >= MATRIX.GetNumberOfRows() || jj >= MATRIX.GetNumberOfColumns())
		{
			std::cout << "provide pairs of values that match the lengths of the sequences:" << "\n";
			std::cout << ii << " should be smaller: " << MATRIX.GetNumberOfRows() << "\n";
			std::cout << jj << " should be smaller: " << MATRIX.GetNumberOfColumns() << "\nbye\n\n";
			exit( -1);
		}
		//std::cout << "anchors: " << ii << "  " << jj << "  " << shift_penalty << "\n";

		for( size_t i = 0; i < MATRIX.GetNumberOfRows(); ++i)
		{
			for( size_t j = 0; j < MATRIX.GetNumberOfColumns(); ++j)
			{
				if
				(
				    ( i  < ii && j  > jj)
				      || ( i  > ii && j  < jj)
				      ||  ( i == ii && j  < jj)
				      || ( i  < ii && j == jj)
				)
				{
					MATRIX( i, j).AddValue( shift_penalty);
				//	std::cout << shift_penalty << "  ";
				}
//				else
//				{
//					std::cout << "0  "; debug
//				}
			}
			//std::cout << "\n"; debug
		}
	}
	return anchors;

	//Kamil was playing here
//	for( size_t i = 0; i < MATRIX.GetNumberOfRows(); ++i)
//	{
//		for( size_t j = 0; j < MATRIX.GetNumberOfColumns(); ++j)
//		{
//				std::cout << MATRIX( i, j).GetValue() << " ";
//			//	std::cout << shift_penalty << "  ";
//		}
//		std::cout << "\n";
//	}
//	//Kamil finished playing here
}

bool
Contains( const std::vector< Triplet< int, int, double> > &ANCHORS, int A, int B)
{
	bool
		contains = false;
	for( std::vector< Triplet< int, int, double> >::const_iterator itr = ANCHORS.begin(); itr != ANCHORS.end(); ++itr )
	{
		if( itr->First() == A && itr->Second() == B){ return true;}
	}
	return contains;
}

bool
ContainsFirst( const std::vector< Triplet< int, int, double> > &ANCHORS, int A)
{
	bool
		contains = false;
	for( std::vector< Triplet< int, int, double> >::const_iterator itr = ANCHORS.begin(); itr != ANCHORS.end(); ++itr )
	{
		if( itr->First() == A){ return true;}
	}
	return contains;
}

bool
ContainsSecond( const std::vector< Triplet< int, int, double> > &ANCHORS, int B)
{
	bool
		contains = false;
	for( std::vector< Triplet< int, int, double> >::const_iterator itr = ANCHORS.begin(); itr != ANCHORS.end(); ++itr )
	{
		if( itr->Second() == B){ return true;}
	}
	return contains;
}



#endif /* ANCHOR_H_ */
