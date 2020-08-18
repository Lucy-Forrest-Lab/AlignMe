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
//!  The ScoreProfileSimilarity returns a simple sum over the differences
//! in a profile.
//! NOTE : NOT USED, NOT IMPLEMENTED!
//! TODO: CLEANUP: there are functions in here that should be elsewhere!
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
  size_t    m_IdInProfile;

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
    DebugWrite( __FUNCTION__ << "  " << d);
    return d;
  }

  virtual std::ostream &Write( std::ostream &STREAM) const
    {
      STREAM << "ScoreProfileSimilarity::Write()" << "\n";
      STREAM << "profile-id: " << m_IdInProfile << "\n";
      return STREAM;
    }
};



// VERSION 1.1 todo: change location!
double PdsScore( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND, const double &GAP_EXTENSION_PENALTY)
{
	double score( 0.0);
	int aligned_count( 0);
	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first == std::numeric_limits< int>::max() || itr->second == std::numeric_limits< int>::max())
		{
			score += pow( GAP_EXTENSION_PENALTY, 2.0); //TODO: use this??? square this?
		}
		else
		{
			std::vector< double>
				first_profile( FIRST[ itr->first].GetProfiles()),
				second_profile( SECOND[ itr->second].GetProfiles());
			for( std::vector< double>::const_iterator first_itr = first_profile.begin(), second_itr = second_profile.begin(); first_itr != first_profile.end() && second_itr != second_profile.end(); ++first_itr, ++second_itr)
			{
				score += /* weight * */ pow( *first_itr - *second_itr, 2.0); // TODO: introduce weights!
			}
			++aligned_count;
		}
	}
	if( aligned_count != 0)
	{
		return sqrt( score / aligned_count);
	}
	return std::numeric_limits< double>::max();
}


// VERSION 1.1 todo: change location!
double PearsonCorrelation( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND, const double &GAP_EXTENSION_PENALTY, const size_t ID = 0, const bool IS_GAPPED = false)
{
	double
		sum_x( 0.0),
		sum_y( 0.0),
		sum_x_square( 0.0),
		sum_y_square( 0.0),
		sum_x_times_y( 0.0),
		length( ALIGNMENT.size());

	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( ( itr->first == std::numeric_limits< int>::max() || itr->second == std::numeric_limits< int>::max()) && IS_GAPPED)
		{
		}
		else
		{
			double
				x( FIRST[ itr->first].GetProfiles()[ ID]),
				y( SECOND[ itr->second].GetProfiles()[ ID]);

			sum_x += x;
			sum_y += y;
			sum_x_square += pow( x, 2.0);
			sum_y_square += pow( y, 2.0);
			sum_x_times_y += x * y;
		}
	}

	return ( length * sum_x_times_y - ( sum_x * sum_y))
		/ ( sqrt( length * sum_x_square - pow( sum_x, 2.0)) * sqrt( length * sum_y_square - pow( sum_y, 2)));
}


#endif
