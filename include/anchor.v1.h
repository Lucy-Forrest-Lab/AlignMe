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
#include "triplet.t.h"

void
BuildAnchorMatrix( std::vector< Triplet< int, int, double> > &LIST, Matrix< DynamicProgrammingMatrixElement> &MATRIX, const int SHIFT = 1)
{
//	std::cout << __FUNCTION__ << ": BE AWARE: IDs are reduced by one, pass them starting with ID = 1 not 0" << std::endl;
	double
		shift_penalty;
	size_t
		ii, jj;

	for( size_t i = 0; i < LIST.size(); ++i)
	{
		ii = LIST[i].First();
		jj = LIST[i].Second();
		shift_penalty = LIST[i].Third();
		ii -= SHIFT;
		jj -= SHIFT;

		if( ii >= MATRIX.GetNumberOfRows() || jj >= MATRIX.GetNumberOfColumns())
		{
			std::cout << "provide pairs of values that match the lengths of the sequences:" << "\n";
			std::cout << ii << " should be smaller: " << MATRIX.GetNumberOfRows() << "\n";
			std::cout << jj << " should be smaller: " << MATRIX.GetNumberOfColumns() << "\n";
			exit( -1);
		}
		//std::cout << "anchors: " << ii << "  " << jj << "  " << shift_penalty << "\n";

		for( size_t i = 0; i < MATRIX.GetNumberOfRows(); ++i)
		{
			for( size_t j = 0; j < MATRIX.GetNumberOfColumns(); ++j)
			{
				if
				(
						   ( i  < ii && j  >= jj)
						|| ( i  >= ii && j  < jj)
//						|| ( i == ii && j  < jj)
//						|| ( i  < ii && j == jj)
				)
				{
					MATRIX( i, j).AddValue( shift_penalty);
				//	std::cout << shift_penalty << "  ";
				}
				else
				{
				//	std::cout << "0  "; debug
				}
			}
			//std::cout << "\n"; debug
		}
	}

	//Kamil was playing here
	for( size_t i = 0; i < MATRIX.GetNumberOfRows(); ++i)
	{
		for( size_t j = 0; j < MATRIX.GetNumberOfColumns(); ++j)
		{
				std::cout << MATRIX( i, j).GetValue() << " ";
			//	std::cout << shift_penalty << "  ";
		}
		std::cout << "\n";
	}
	//Kamil finished playing here
}



#endif /* ANCHOR_H_ */
