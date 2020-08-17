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
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef SCORE_SEQUENCE_SIMILARITY_H
#define SCORE_SEQUENCE_SIMILARITY_H

#include <set>

#include "function.t.h"
#include "amino_acid.h"
#include "macro_functions_read_write.h"


class ScoreSequenceSimilarity
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
 protected:
  // for a pair of aminoacids return substitution score
  std::map< std::string, double>     m_AAPairToScore;

 public:
  //! default constructor
  ScoreSequenceSimilarity()
    :  m_AAPairToScore()
    {}


    //! construct from file name of similarity matrix
    ScoreSequenceSimilarity( const std::map< std::string, double> &MAP)
      : m_AAPairToScore( MAP)
    {
		DebugWrite( __FUNCTION__);
    }

    ScoreSequenceSimilarity( const ScoreSequenceSimilarity &ORIG)
    : m_AAPairToScore( ORIG.m_AAPairToScore)
      {
		DebugWrite( __FUNCTION__ << " copy constructor");
      }

    //! virtual destructor
    virtual ~ScoreSequenceSimilarity()
      {}
    
    //
    virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA);

    const std::map< std::string, double> &GetSubstitutionMap() const
    {
    	return m_AAPairToScore;
    }

    void SetSubstitutionMap( const std::map< std::string, double> &MAP)
    {
    	m_AAPairToScore = MAP;
    }



	virtual int GetClassID() const
	{
		return e_SequenceSimilarity;
	}

  
  virtual std::ostream &Write( std::ostream &STREAM) const
    {
      STREAM << "ScoreSequenceSimilarity::Write()" << "\n";
#ifdef DEBUG
      STREAM << m_AAPairToScore;
#endif
      return STREAM;
    }

};



double MatrixValueRange( const std::map< std::string, double> &MATRIX);


std::map< std::string, double>
ScaleSubstitutionMatrixToUnity( const std::map< std::string, double> &MATRIX, const double &MIN_CUTOFF = std::numeric_limits< double>::max());

// read the substitution matrix into the map
std::map< std::string, double>
ReadSubstitutionMatrix( const std::string &FILE, const std::set< char> &DEFINED_AMINO_ACIDS);


#endif


