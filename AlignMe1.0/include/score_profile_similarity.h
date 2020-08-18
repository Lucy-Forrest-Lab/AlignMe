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
#include "amino_acid.h"

#include <list>


class ScoreProfileSimilarity
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
 protected: 
  size_t    m_IdInProfile;    //!< ID of the profile in the GeneralizedAminoAcid that an object of this class refers to, needed in operator()

 public:
  
  ScoreProfileSimilarity()
    : m_IdInProfile( std::numeric_limits< size_t>::max())
    {}

  ScoreProfileSimilarity( const size_t &ID)
    : m_IdInProfile( ID)
    {
      DebugWrite( __FUNCTION__);
    }

  virtual ~ScoreProfileSimilarity(){}


  virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
  {
    double d = -fabs( AA.first.GetProfiles()[ m_IdInProfile] - AA.second.GetProfiles()[ m_IdInProfile]);
 //   DebugWrite( __FUNCTION__ << "  " << d);
    return d;
  }

  virtual std::ostream &Write( std::ostream &STREAM) const
    {
      STREAM << "ScoreProfileSimilarity::Write()" << "\n";
      STREAM << "profile-id: " << m_IdInProfile << "\n";
      return STREAM;
    }
};


int NrAligned( const std::list< std::pair< int, int> > &ALIGNMENT)
{
	int count = 0;
	for( std::list< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max())
		{
			++count;
		}
	}
	return count;
}


int LengthAligned( const std::list< std::pair< int, int> > &ALIGNMENT)
{
	int length = 0;
	for( std::list< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		++length;
	}
	return length;
}


#endif
