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

#include <list>
#include <iterator>

#include "function.t.h"
#include "amino_acid.h"
#include "string_functions.h"


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
    ScoreSequenceSimilarity( const std::string &FILENAME, const std::set< char> &DEFINED_AMINO_ACIDS)
      : m_AAPairToScore()
    {
		DebugWrite( __FUNCTION__);
		Read( FILENAME, DEFINED_AMINO_ACIDS);
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
    virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
  {
    	DebugWrite( __FUNCTION__);
    std::string 
      amino_acid_pair( std::string( 1, AA.first.GetType()) + std::string( 1, AA.second.GetType()));
    return m_AAPairToScore[ amino_acid_pair];  // ask substitution matrix/map for value of AA pair 'AB'
  }

    // read the substitution matrix into the map
    void Read( const std::string &FILE, const std::set< char> &DEFINED_AMINO_ACIDS)
    {
    	std::ifstream read( FILE.c_str());
    	if( !read)
    	{
    		std::cerr << "ERROR: File <" << FILE << "> could not be opened! " << "\n";
    		exit( -1);
    	}
		  DebugWrite( __FUNCTION__);

		  char first_aa;

		  double value;

		  std::string
			  str,
			  line;

		  std::getline( read, line);

		  DebugWrite( "read line: " << line);


		  std::vector< std::string> amino_acids  = mystr::SplitString( line);

#ifdef DEBUG
		  std::copy( amino_acids.begin(), amino_acids.end(), std::ostream_iterator< std::string>( std::cout, "\n"));
#endif

#ifdef SECURE
		  for( std::set< char>::const_iterator itr = DEFINED_AMINO_ACIDS.begin(); itr != DEFINED_AMINO_ACIDS.end(); ++itr)
		  {
			  if( std::find( amino_acids.begin(), amino_acids.end(), std::string( 1, *itr)) == amino_acids.end())
			  {
				  std::cerr << "ERROR: the amino acid <" << *itr << "> in the sequence input is not found in the substitution matrix " << FILE << "\n";
				  exit( -1);
			  }
		  }
#endif

		  int nr_amino_acids = amino_acids.size();


		  for( int i = 0; i < nr_amino_acids; ++i) // lines
		  {


#ifdef SECURE
			  str.clear();

			  read >> str;

			  if( str.size() == 0)
			  {
				  continue;
			  }
			  else if( str.size() > 1)
			  {
				  std::cerr << "ERROR: Amino acid identifier only in one letter code!: <" << str << ">" << "\n";
				  exit( -1);
			  }

			  first_aa = str[0];

			  if( first_aa != amino_acids[i][0])
			  {
				  std::cerr << "ERROR: The substitution matrix "<< FILE <<" is not symmetric "<< "\n";
				  std::cerr << "ERROR: Only amino acid codes are allowed in the first row and the first column of the matrix-file" << "\n";
				  std::cerr << "ERROR: At position "<< i+1 << " the program found <" << amino_acids[i][0] << "> in the first row but <" << first_aa << "> in the first column" << "\n";
				  exit( -1);
			  }

#else
			  read >> first_aa;
#endif


			  for( int j = 0; j < nr_amino_acids; ++j)
			  {
#ifdef SECURE
				  str.clear();
				  read >> str;
				  if( !mystr::IsNumerical( str))
				  {
					  std::cerr << "ERROR: Expected a numerical value in the substitution matrix "<< FILE <<" but found: <" << str << "> at position ["<< i+1 <<"]["<<j+1<<"]"  << "\n";
					  exit( -1);
				  }
				  value = mystr::ConvertStringToNumericalValue< double>( str);
#else
				  read >> value;
#endif

				  std::string amino_acid_pair;
				  amino_acid_pair += first_aa;
				  amino_acid_pair += amino_acids[ j];
				  m_AAPairToScore[ amino_acid_pair] = value;

			  }
		  }

		  read.close();
		  read.clear();
    }


  
  virtual std::ostream &Write( std::ostream &STREAM) const
    {
      STREAM << "ScoreSequenceSimilarity::Write()" << "\n";
#ifdef DEBUG
      STREAM << m_AAPairToScore;
#endif
      return STREAM;
    }


  double MatrixValueRange()
  {
	  double
		  min( std::numeric_limits< double>::max()),
		  max( -1.0 * min);

	  for( std::map< std::string, double>::const_iterator itr = m_AAPairToScore.begin(); itr != m_AAPairToScore.end(); ++itr)
	  {
		  if( itr->second > max)
		  {
			  max = itr->second;
		  }
		  if( itr->second < min)
		  {
			  min = itr->second;
		  }
	  }
	  return max - min;
  }


};

#endif


