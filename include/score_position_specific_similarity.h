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


#ifndef SCORE_POSITION_SPECIFIC_SIMILARITY_H
#define SCORE_POSITION_SPECIFIC_SIMILARITY_H

#include <list>
#include <set>
#include <iterator>

#include "function.t.h"
#include "amino_acid.h"
#include "string_functions.h"
#include "macro_functions_read_write.h"


class ScorePositionSpecificSimilarity
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
public:
	//! default constructor
	ScorePositionSpecificSimilarity()
	{}


	ScorePositionSpecificSimilarity( const ScorePositionSpecificSimilarity &ORIG)
    {
    	DebugWrite( __FUNCTION__ << " copy constructor");
    }

    //! virtual destructor
    virtual ~ScorePositionSpecificSimilarity()
    {}
    
    //
    virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
    {
    	DebugWrite( __FUNCTION__);


    	return 0.5 *( AA.first.GetPSSMValue(AA.second.GetType()) +  AA.second.GetPSSMValue(AA.first.GetType()) );
//      	return 0.5 * (
//      			AA.first.GetSeqPtr()->GetPositionSpecificSimilarityMatrix()[ AA.first.GetSeqID()][ AA.second.GetType()]
//      		    + AA.second.GetSeqPtr()->GetPositionSpecificSimilarityMatrix()[ AA.second.GetSeqID()][ AA.first.GetType()]
 //   	return 0.5 * ( AA.first.GetPSSScore( AA.second.GetType()) + AA.second.GetPSSScore( AA.first.GetType()));
    }



//#ifndef SEQPTR
//    // read the substitution matrix into the map
//    std::vector< std::map< char, int> > Read( const std::string &FILE, const std::set< char> &DEFINED_AMINO_ACIDS) //, const Sequence *SEQ)
//    {
//
//    	std::string
////			defined = "ARNDCQEGHILKMFPSTWYVBZX",
//			defined = "XAXCDEFGHIKLMNPQRSTVWXYXXX",
//			sequence,
//			line,
//			tmp;
//
//    	std::cout << __FUNCTION__ ;
//    	std::cout << "| defined aa's: " << defined << " (" << defined.size() << ")" << std::endl;
//
////    	std::cout << "min short " << std::numeric_limits< short int>::min() << std::endl;exit( 1);
//    	std::ifstream
//			read( FILE.c_str());
//
//    	if( !read)
//    	{
//    		std::cerr << "ERROR: File <" << FILE << "> could not be opened! " << "\n";
//    		exit( -1);
//    	}
//
//
//    	int
//			min = std::numeric_limits< int>::max(),
//			max = -std::numeric_limits< int>::max(),
//			nr_aa;
//
//    	read >> nr_aa;
//
//    	read >> sequence;
//
//    	std::vector< std::map< char, int> >
//			scores( nr_aa);
//
//    	std::cout << "sequence: " << sequence << std::endl;
//    	std::cout << "size: "<< nr_aa << std::endl;
//
//    	if( nr_aa != (int) sequence.size())
//    	{
//			std::cout << "sizes do not match: " << nr_aa << " != size(" << sequence << ")" << std::endl;
//			exit( -1);
//    	}
//
//    	read >> tmp >> tmp >> tmp >> tmp;
//    	read >> tmp >> tmp >> tmp >> tmp;
//    	read >> tmp >> tmp >> tmp >> tmp;
//    	std::getline( read, line); // get rid of newline
//
//    	std::string::const_iterator
//			seq_itr = sequence.begin(),
//			def_itr;
//
//    	std::vector< int>::const_iterator
//			val_itr;
//
//    	std::vector< std::map< char, int> >::iterator
//			score_itr = scores.begin();
//
//    	for( int i = 0; i < nr_aa; ++i, ++seq_itr, ++score_itr)
//    	{
//    		line.clear();
//
//    		std::getline( read, line);
//
//    		DebugWrite( "read line: " << line);
////    		std::cout << "read line "<<i<<": " << line << std::endl;
//
//    		std::vector< int>
//				values = mystr::SplitAndConvertString< int>( line);
//
//#ifdef SECURE
//    		if( (int) values.size() < (int) defined.size())
//    		{
//    			std::cout << "too few values given, expected at least " << defined.size() << " got: " << values << std::endl;
//    			exit( -1);
//    		}
//#endif
//
//    		def_itr = defined.begin();
//    		val_itr = values.begin();
//			for( int j = 0; j < (int) defined.size(); ++j, ++def_itr, ++val_itr)
//			{
////				std::cout << *val_itr << " ";
//				scores[ i][ *def_itr] = *val_itr;
//				if( min > *val_itr && *val_itr > std::numeric_limits< short int>::min())
//				{
//					min = *val_itr;
//				}
//				if( max < *val_itr)
//				{
//					max = *val_itr;
//				}
//			}
////			std::cout << std::endl;
//
//    	}
//    	std::cout << "range of values in PositionSpecificSimilarity matrix: " << min << " : " << max << std::endl;
//    	read.close();
//    	read.clear();
//    	return scores;
//    }
//#endif


	virtual int GetClassID() const
	{
		return e_PositionSpecificSimilarity;
	}

  
  virtual std::ostream &Write( std::ostream &STREAM) const
    {
      STREAM << "ScorePositionSpecificSimilarity::Write()" << "\n";
//#ifdef DEBUG
//      STREAM << m_PSIScore;
//#endif
      return STREAM;
    }

};





class MinScorePositionSpecificSimilarity
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
public:
	//! default constructor
	MinScorePositionSpecificSimilarity()
	{}


	MinScorePositionSpecificSimilarity( const MinScorePositionSpecificSimilarity &ORIG)
    {
    	DebugWrite( __FUNCTION__ << " copy constructor");
    }

    //! virtual destructor
    virtual ~MinScorePositionSpecificSimilarity()
    {}

    //
    virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
    {
    	DebugWrite( __FUNCTION__);


    	return std::min( AA.first.GetPSSMValue(AA.second.GetType()), AA.second.GetPSSMValue(AA.first.GetType()) );
   }

	virtual int GetClassID() const
	{
		return e_PositionSpecificSimilarity;
	}


  virtual std::ostream &Write( std::ostream &STREAM) const
    {
      STREAM << "MinScorePositionSpecificSimilarity::Write()" << "\n";
      return STREAM;
    }

};



class ProfileScorePositionSpecificSimilarity
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
public:
	//! default constructor
	ProfileScorePositionSpecificSimilarity()
	{}


	ProfileScorePositionSpecificSimilarity( const ProfileScorePositionSpecificSimilarity &ORIG)
    {
    	DebugWrite( __FUNCTION__ << " copy constructor");
    }

    //! virtual destructor
    virtual ~ProfileScorePositionSpecificSimilarity()
    {}

    //
    virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
    {
    	DebugWrite( __FUNCTION__);

       double d = 0;

        for (std::map< char, int>::const_iterator itr=AA.first.GetPSSM().begin();  itr != AA.first.GetPSSM().end() ; ++itr )
        {
        	d -=  double( abs( itr->second - AA.second.GetPSSMValue(itr->first) ) );
        }

        return d / double ( AA.first.GetPSSM().size());
   }

	virtual int GetClassID() const
	{
		return e_PositionSpecificSimilarity;
	}


  virtual std::ostream &Write( std::ostream &STREAM) const
    {
      STREAM << "ProfileScorePositionSpecificSimilarity::Write()" << "\n";
      return STREAM;
    }

};



class fdotfScorePositionSpecificSimilarity
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
public:
	//! default constructor
	fdotfScorePositionSpecificSimilarity()
	{}


	fdotfScorePositionSpecificSimilarity( const fdotfScorePositionSpecificSimilarity &ORIG)
    {
    	DebugWrite( __FUNCTION__ << " copy constructor");
    }

    //! virtual destructor
    virtual ~fdotfScorePositionSpecificSimilarity()
    {}

    //
    virtual double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
    {
    	DebugWrite( __FUNCTION__);

       double d = -1;

        for (std::map< char, int>::const_iterator itr=AA.first.GetPSSM().begin();  itr != AA.first.GetPSSM().end() ; ++itr )
        {
        	d +=  double( itr->second * AA.second.GetPSSMValue(itr->first));
        }

        return 0.0001 * d ;
   }

	virtual int GetClassID() const
	{
		return e_PositionSpecificSimilarity;
	}


  virtual std::ostream &Write( std::ostream &STREAM) const
    {
      STREAM << "fdotfScorePositionSpecificSimilarity::Write()" << "\n";
      return STREAM;
    }

};



#endif


