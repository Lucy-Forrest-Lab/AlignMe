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

#include "definitions.h"

#include "matrix.t.h"
#include "dynamic_programing_matrix_element.h"
#include "triplet.t.h"

std::vector< Triplet< int, int, double> >
BuildAnchorMatrix( std::istream &STREAM, Matrix< DynamicProgrammingMatrixElement> &MATRIX, const double &PENALTY = -1.0);


bool
Contains( const std::vector< Triplet< int, int, double> > &ANCHORS, int A, int B);


bool
ContainsFirst( const std::vector< Triplet< int, int, double> > &ANCHORS, int A);


bool
ContainsSecond( const std::vector< Triplet< int, int, double> > &ANCHORS, int B);

#endif /* ANCHOR_H_ */
