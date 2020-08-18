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


#ifndef SCORE_SEQUENCE_SIMILARITY_PROFILE_DEPENDENT_H
#define SCORE_SEQUENCE_SIMILARITY_PROFILE_DEPENDENT_H

#include "function.t.h"
#include "amino_acid.h"
#include "macro_functions_read_write.h"


class ScoreProfileDependentSequenceSimilarity
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
 protected:
  // for a pair of aminoacids return substitution score
	 std::map< std::string, double>     m_Upper;
	 std::map< std::string, double>     m_Lower;
	 size_t						        m_ProfileID;
	 double							    m_Threshold;

 public:
  //! default constructor
  ScoreProfileDependentSequenceSimilarity()
    :  m_Upper(), m_Lower(), m_ProfileID(),m_Threshold()
    {}


    //! construct from file name of similarity matrix
    ScoreProfileDependentSequenceSimilarity( const std::map< std::string, double> &UPPER, const std::map< std::string, double> &LOWER, const int &PROFILE_ID, const double &THRESHOLD)
      : m_Upper( UPPER),
        m_Lower( LOWER),
        m_ProfileID( PROFILE_ID),
        m_Threshold( THRESHOLD)
    {
		DebugWrite( __FUNCTION__);
    }

    ScoreProfileDependentSequenceSimilarity( const ScoreProfileDependentSequenceSimilarity &ORIG)
    : m_Upper( ORIG.m_Upper),
      m_Lower( ORIG.m_Lower),
      m_ProfileID( ORIG.m_ProfileID),
      m_Threshold( ORIG.m_Threshold)
      {
		DebugWrite( __FUNCTION__ << " copy constructor");
      }

    //! virtual destructor
    virtual ~ScoreProfileDependentSequenceSimilarity()
    {}
    
    //
    virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA);


	virtual int GetClassID() const;

  
	virtual std::ostream &Write( std::ostream &STREAM) const;



};


#endif


