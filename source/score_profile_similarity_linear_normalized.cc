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
//! @date: 27.5.2011
/////////////////////////////////////////////////////////////////////////

#include "../include/score_profile_similarity_linear_normalized.h"

#include <cmath>
#include "../include/macro_functions_read_write.h"



ScoreProfileSimilarityLinearNormalized::ScoreProfileSimilarityLinearNormalized()
: m_IdInProfile( std::numeric_limits< size_t>::max()),
  m_DeltaMaxInv( 0)
{}

ScoreProfileSimilarityLinearNormalized::ScoreProfileSimilarityLinearNormalized( const size_t &ID, const double &DELTA_MAX)
: m_IdInProfile( ID),
  m_DeltaMaxInv( 1.0 / DELTA_MAX)
{
  DebugWrite( __FUNCTION__);
}

ScoreProfileSimilarityLinearNormalized::~ScoreProfileSimilarityLinearNormalized(){}


double ScoreProfileSimilarityLinearNormalized::operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
{
//	std::cout << "id: " << m_IdInProfile << std::endl;
//	std::cout << AA.first.GetProfiles()[ m_IdInProfile] << std::endl;
//	std::cout << AA.second.GetProfiles().size() << std::endl;
	double d = 0 - (m_DeltaMaxInv * fabs( AA.first.GetProfiles()[ m_IdInProfile] - AA.second.GetProfiles()[ m_IdInProfile]));
	//std::cout << d << "_";
	//   DebugWrite( __FUNCTION__ << "  " << d);
	return d;
}

std::ostream &ScoreProfileSimilarityLinearNormalized::Write( std::ostream &STREAM) const
{
  STREAM << "ScoreProfileSimilarityLinearNormalized::Write()" << "\n";
  STREAM << "profile-id: " << m_IdInProfile << "\n";
  STREAM << "inverse-delta-max: " << m_DeltaMaxInv << "\n";
  return STREAM;
}

