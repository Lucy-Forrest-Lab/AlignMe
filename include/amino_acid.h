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
//!  amino_acid.h
//! The GeneralizedAminoAcid class stores all profiles/scales for a
//! residue. It can have an amino acid identity or not.
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef GENERALIZED_AMINO_ACID_H
#define GENERALIZED_AMINO_ACID_H

#include "std_functions.h"
#include "definitions.h"

// A GeneralizedAminoAcid always consists of a m_Type and a m_Profiles

//class Sequence;

class GeneralizedAminoAcid
{
protected:

	char                   m_Type;
	std::vector< double>   m_Profiles;
	unsigned short         m_SeqID;
#ifdef POSITION_SPECIFIC_SIMILARITY
	std::map< char, int>   m_PSSM; //!<  position specific substitution matrix
#endif

public:

	GeneralizedAminoAcid()
	: m_Type('X'),
	m_Profiles(),
	m_SeqID()
#ifdef POSITION_SPECIFIC_SIMILARITY
	,m_PSSM()
#endif
	{}


	// If a certain TYPE is given, than it will be added to the m_Type of the GAA (f.e. type AminoAcid A, G, I and so on)
	GeneralizedAminoAcid( const char &TYPE, const size_t &SIZE = 0)
	: m_Type( TYPE),
	m_Profiles( SIZE, 0.0),
	m_SeqID()
#ifdef POSITION_SPECIFIC_SIMILARITY
	,m_PSSM()
#endif
	{}

	GeneralizedAminoAcid( const GeneralizedAminoAcid &GAA)
	: m_Type( GAA.m_Type),
    m_Profiles( GAA.m_Profiles),
	m_SeqID( GAA.m_SeqID)
#ifdef POSITION_SPECIFIC_SIMILARITY
	,m_PSSM( GAA.m_PSSM)
#endif
{ /*std::cout << __PRETTY_FUNCTION__ << "\n";*/}

	virtual ~GeneralizedAminoAcid(){}

	// returns the m_Type of a certain GAA
	virtual const char &GetType() const
	{
		return m_Type;
	}

	void SetSeqID( const unsigned short &ID)
	{
		m_SeqID = ID;
	}

	unsigned short GetSeqID() const
	{
		return m_SeqID;
	}

#ifdef POSITION_SPECIFIC_SIMILARITY

	virtual void SetPSSM(const char &AA_TYPE, const int &VALUE)
	{
		m_PSSM[AA_TYPE] = VALUE;
	}

	virtual const int &GetPSSMValue(const char &AA_TYPE) const
	{
		return m_PSSM.find(AA_TYPE)->second;  // TODO: UNSAFE: if AA_TYPE is not found the function will crash!
	}

	virtual const std::map< char, int> &GetPSSM() const
	{
		return m_PSSM;
	}


#endif


  // returns whole profile of a certain GAA. There could be several values stores in m_Profile
  // because it is a vector of doubles
  virtual const std::vector< double> &GetProfiles() const
  {
    return m_Profiles;
  }

  // adds a new value to the m_Profile of a certain GAA
  virtual void AddNewProfile( const double &VALUE)
  {
    m_Profiles.push_back( VALUE);
  }

  virtual void SumToProfiles( const std::vector< double> &VALUES)
  {
	  std::vector< double>::iterator prof_itr = m_Profiles.begin();
	  for( std::vector< double>::const_iterator val_itr = VALUES.begin(); val_itr != VALUES.end(); ++val_itr, ++prof_itr)
	  {
		  *prof_itr += *val_itr;
	  }
  }

  virtual void Normalize( const double &FACTOR)
  {
	  for( std::vector< double>::iterator itr = m_Profiles.begin(); itr != m_Profiles.end(); ++itr)
	  {
		  *itr /= FACTOR;
	  }
  }

  std::ostream &Write( std::ostream &STREAM) const
  {
      STREAM << m_Type << "\n";
      STREAM << m_Profiles;
      return STREAM;
  }
  
};


inline
std::ostream &operator << ( std::ostream &STREAM, const GeneralizedAminoAcid &AA)
{
  return AA.Write( STREAM);
}




#endif
