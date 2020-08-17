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
//!  The classical sequence similarity score based on substitution matrices.
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: @date: 27.5.2011
/////////////////////////////////////////////////////////////////////////


#include "../include/score_sequence_similarity_profile_dependent.h"

#include "../include/string_functions.h"

    
//
double ScoreProfileDependentSequenceSimilarity::operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
{
	DebugWrite( __FUNCTION__);
	std::string
	amino_acid_pair( std::string( 1, AA.first.GetType()) + std::string( 1, AA.second.GetType()));
#ifdef SECURE
	if (AA.first.GetProfiles().size() == 0)
	{
		std::cerr << "ERROR: No profile provided! If you want to use the type: 'NormalizedProfileDependentSubstitutionMatrix' then you also have to provide at least one profile (scale or prediction) \n";
		exit(-1);
	}
	if (m_ProfileID >= AA.first.GetProfiles().size())
	{
		std::cerr << "ERROR: Profile number [" << m_ProfileID << "] could not be accessed, because there are only " <<  AA.first.GetProfiles().size() << " profiles in your alignment \n";
		exit(-1);
	}
#endif
	if( AA.first.GetProfiles()[ m_ProfileID] > m_Threshold && AA.second.GetProfiles()[ m_ProfileID] > m_Threshold)
	{
		return m_Upper[ amino_acid_pair];
	}
	return m_Lower[ amino_acid_pair];  // ask substitution matrix/map for value of AA pair 'AB'
}


int ScoreProfileDependentSequenceSimilarity::GetClassID() const
{
	return e_ProfileDependentSequenceSimilarity;
}


std::ostream &ScoreProfileDependentSequenceSimilarity::Write( std::ostream &STREAM) const
{
  STREAM << "ScoreProfileDependentSequenceSimilarity::Write()" << "\n";
#ifdef DEBUG
  STREAM << "upper: \n" << m_Upper << "\nlower:\n" << m_Lower << std::endl;
#endif
  STREAM << "profile-id: " << m_ProfileID << std::endl;
  STREAM << "threshold: " << m_Threshold << std::endl;
  return STREAM;
}


