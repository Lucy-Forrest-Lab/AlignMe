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
//! The ScoreProfileSimilarityLinearNormalized class returns a sum over the differences.
//!
//! TODO: seperate functions from class
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef SCORE_PROFILE_SIMILARITY_LINEAR_NORMALIZED_H
#define SCORE_PROFILE_SIMILARITY_LINEAR_NORMALIZED_H

#include "function.t.h"
#include "sequence.h"

#include <list>


class ScoreProfileSimilarityLinearNormalized
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
protected:
	size_t    m_IdInProfile;    //!< ID of the profile in the GeneralizedAminoAcid that an object of this class refers to, needed in operator()
	double    m_DeltaMaxInv;

public:
  
	ScoreProfileSimilarityLinearNormalized();

	ScoreProfileSimilarityLinearNormalized( const size_t &ID, const double &DELTA_MAX);

	virtual ~ScoreProfileSimilarityLinearNormalized();

	virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA);

	virtual std::ostream &Write( std::ostream &STREAM) const;

	virtual int GetClassID() const
	{
		return e_LinearNormalizedProfileSimilarity;
	}
}; // end class


#endif
