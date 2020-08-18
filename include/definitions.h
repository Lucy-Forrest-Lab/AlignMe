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

// defining SECURE will enable a lot of security checks:
#define SECURE

// lots and lots of very useful output:
//#define DEBUG


// basic definition of an amino acid sequence as vector of generalized amino acids:
#include <vector>
#include "amino_acid.h"
typedef std::vector< GeneralizedAminoAcid>  AASequence;


#endif // DEFINITIONS_H
