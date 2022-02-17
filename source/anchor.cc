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

#include "anchor.h"


std::vector< Triplet< int, int, double> >
BuildAnchorMatrix( std::istream &STREAM, Matrix< DynamicProgrammingMatrixElement> &MATRIX, const double &PENALTY)
{
	// anchors need to be given in usual vector notation
	double
		shift_penalty;
	size_t
		nr, anchor_pos_1, anchor_pos_2;
	std::vector< Triplet<int,int,double> >
		anchors;
	STREAM >> nr;
	for( size_t r = 0; r < nr; ++r)
	{
		STREAM >> anchor_pos_1 >> anchor_pos_2 >> shift_penalty;
		if( anchor_pos_1 == 0 || anchor_pos_2 == 0)
		{
			std::cout << "Counting of the amino acids start with '1'." << std::endl;
			std::cout << "anchors for the first residue have to be indicated by '1' " << std::endl;
			exit( 2);
		}
		anchors.push_back( Triplet<int,int,double>( anchor_pos_1, anchor_pos_2, shift_penalty));
		shift_penalty *= -1;  // using positive values on website, here they need to be negative, similar as a potential

		if( anchor_pos_1 >= MATRIX.GetNumberOfRows() || anchor_pos_2 >= MATRIX.GetNumberOfColumns())
		{
			std::cout << "provide pairs of values that match the lengths of the sequences:" << "\n";
			std::cout << anchor_pos_1 << " should not exceed the max value: " << MATRIX.GetNumberOfRows() - 1 << "\n";
			std::cout << anchor_pos_2 << " should not exceed the max value: " << MATRIX.GetNumberOfColumns() - 1 << "\nbye\n\n";
			exit( 2);
		}

		for( size_t i = 0; i < MATRIX.GetNumberOfRows(); ++i)
		{
			for( size_t j = 0; j < MATRIX.GetNumberOfColumns(); ++j)
			{
				if
				(
				    ( i  < anchor_pos_1 && j  > anchor_pos_2)
				      || ( i  > anchor_pos_1 && j  < anchor_pos_2)
				      ||  ( i == anchor_pos_1 && j  < anchor_pos_2)
				      || ( i  < anchor_pos_1 && j == anchor_pos_2)
				)
				{
					MATRIX( i, j).AddValue( shift_penalty);
				}
			}
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

