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
//!  A collection of functions to evaluate the similarity of two profiles.
//! The ScoreProfileSimilarity class returns a sum over the differences.
//!
//! TODO: seperate functions from class
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef SCORE_PROFILE_SIMILARITY_H
#define SCORE_PROFILE_SIMILARITY_H

#include "function.t.h"
#include "sequence.h"

#include <list>


class ScoreProfileSimilarity
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
 protected: 
  size_t    m_IdInProfile;    //!< ID of the profile in the GeneralizedAminoAcid that an object of this class refers to, needed in operator()

 public:
  
  ScoreProfileSimilarity();

  ScoreProfileSimilarity( const size_t &ID);

  virtual ~ScoreProfileSimilarity();

  virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA);

  virtual std::ostream &Write( std::ostream &STREAM) const;

}; // end class



int NrAligned( const std::vector< std::pair< int, int> > &ALIGNMENT);

int LengthAligned( const std::vector< std::pair< int, int> > &ALIGNMENT);

// VERSION 1.1================================ till end of file !
//This score is good when aligning sequences of similar length and when profiles match well. Do not care explicitly about gaps.
double PdsScoreAlignmentLength
(
		const std::vector< std::pair< int, int> > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND,
		const double &RATIO = 0.4
);

//This score is good when profiles match well and do not care of differences in length. Do not care explicitly about gaps.
double PdsScoreAligned
(
		const std::vector< std::pair< int, int> > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND,
		const double &RATIO = 0.1
);


//This score is good when profiles match well and few gaps were found.
double PdsScoreAlignedAndGaps
(
		const std::vector< std::pair< int, int> > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND,
		const double &RATIO = 0.1
);


double PdsScoreAlignedMin( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND);

double PdsScoreAlignedMax( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND);

//This score is good when aligning sequences of similar length (because of gaps!) and when profiles match well. Mean value for gaps.
double PdsScore
(
		const std::vector< std::pair< int, int> > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND,
		const double &GAP_EXTENSION_PENALTY = 0.0,
		const double &RATIO = 0.1
);


double PdsFixedGaps
(
		const std::vector< std::pair< int, int> > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND,
		const double &GAP_EXTENSION_PENALTY
);


std::vector< double> CalculatePearsonCorrelations
(
		const std::vector< std::pair< int, int> > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND,
		const double &RATIO = 0.2
);


#endif
