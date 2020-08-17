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
//!  General includes, typedefs and macro definitions are placed here.
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

// defining SECURE will enable a lot of security checks:
#define SECURE

// lots and lots of very useful output:
//#define DEBUG


// use when compiling msa_optimizer.cpp:
#define POSITION_SPECIFIC_SIMILARITY

#define Alignment std::vector< std::pair< int, int> >

#define ShPtr boost::shared_ptr

#define Profile std::vector< std::pair< double, double> >

// for profile values
#define g_GAP 888888.0

// for alignment indices
#define g_GapID std::numeric_limits< int>::max()

#define MSA std::vector< Sequence>

enum SimilarityType { e_SequenceSimilarity, e_ProfileSimilarity, e_ProbabilityProfileSimilarity, e_PositionSpecificSimilarity, e_ProfileDependentSequenceSimilarity, e_LinearNormalizedProfileSimilarity};


#endif // DEFINITIONS_H
