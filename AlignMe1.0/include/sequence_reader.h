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
//!  The main reading functions.
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef SEQUENCE_READER_H
#define SEQUENCE_READER_H

#include "sliding_window.h"
#include "string_functions.h"




//! read profiles and scales into two sets of mul
void  ReadProfilesAndScales
(
		std::vector< Sequence> &SEQUENCES,
		const std::string &SCORE_AND_PROFILE_FILE,
		std::set<char> &DEFINED_AMINO_ACIDS
)
{
	DebugWrite(  __FUNCTION__ );
	std::string
		ending,
		type,
		file,
		tmp,
		sliding_window_type;

	double
		weight;

	int
		window_size,
		column,
		header_lines;

	std::ifstream
		read( SCORE_AND_PROFILE_FILE.c_str());

	if( !read)
	{
		std::cerr << "ERROR: The score_and_profile_file [" << SCORE_AND_PROFILE_FILE << "] was not found" << "\n";
		exit( -1);
	}

	// read the score and profile file
	// for each command-line this loop will be executed
  while( read)
  {
	  type.clear();
	  tmp.clear();
#ifdef SECURE
	  read >> tmp >> tmp;
	  if (tmp == "")
	  {
		  continue;

	  }
	  if( !mystr::IsNumerical( tmp) || (weight = mystr::ConvertStringToNumericalValue<int>( tmp) < 0))
	  {
		  std::cerr << "ERROR: weight has to be positive number or zero, found: <" << tmp << ">" << "\n";
		  exit( -1);
	  }
#else
	  read >> tmp >> weight;
#endif

	  read >> tmp >> type;

		//  weight similaritytype (f.e. BLOSUM or PAM)
	  if( type == "SequenceSimilarity")
	  {
		  read >> tmp >> tmp;
	  }
		// weight similaritytype firstprofile_file secondprofile_file
		// FIRST = first sequence  SECOND = second sequence
	  else if( type == "ProfileSimilarity")
	  {
		  read >> tmp >> ending;
		  for( std::vector< Sequence>::iterator itr = SEQUENCES.begin(); itr != SEQUENCES.end(); ++itr)
		  {
			  file = itr->GetFileName();
			  file = file.substr( 0, file.find( "."));
			  file += ending;
			  ReadProfile( *itr, file);
		  }
	  }
	  else if( type == "UniversalProfileSimilarity")
	  {

#ifdef SECURE
		  read >> tmp >> tmp;
		  if( !mystr::IsNumerical( tmp) || ( column = mystr::ConvertStringToNumericalValue< size_t>( tmp)) < 0)
		  {
			  std::cerr<< "ERROR: column value for UniveralProfileSimilarity has to be a positive number or zero, found: <"<< tmp<< "> in " <<  SCORE_AND_PROFILE_FILE << "\n";
			  exit( -1);
		  }
		  read >> tmp >> tmp;
		  if( !mystr::IsNumerical( tmp) || ( header_lines = mystr::ConvertStringToNumericalValue< size_t>( tmp)) < 0)
		  {
			  std::cerr<< "ERROR: header_lines value for UniveralProfileSimilarity has to be a positive number or zero, found: <"<< tmp<< "> in " <<  SCORE_AND_PROFILE_FILE << "\n";
			  exit( -1);
		  }
#else
		  read >> tmp >> column;
		  read >> tmp >> header_lines;
#endif
		  read >> tmp >> ending;

		  for( std::vector< Sequence>::iterator itr = SEQUENCES.begin(); itr != SEQUENCES.end(); ++itr)
		  {
			  file = itr->GetFileName();
			  file = file.substr( 0, file.find( "."));
			  file += ending;
			  ReadUniversalProfile( *itr, file, column, header_lines);
		  }
	  }
		// weight similaritytype scale_file windowtype windowsize
		// scale is f.e. KD or Eisenberg
		// checks also if there are two sequences (FIRST and SECOND) that can be compared
	  else if( type == "ScaleSimilarity")
	  {
		  tmp.clear();
		  read >> tmp >> file;
#ifdef SECURE
		  if (file == "")
		  {
			  std::cerr << "ERROR: You have to provide a scale if you want to do an alignment based on scale similarity"<<"\n";
			  exit(-1);
		  }
#endif
		  tmp.clear();
		  read >> tmp >> sliding_window_type;
		  read >> tmp >> window_size;

		  for( std::vector< Sequence>::iterator itr = SEQUENCES.begin(); itr != SEQUENCES.end(); ++itr)
		  {
			  ReadScale( *itr, file, sliding_window_type, window_size, DEFINED_AMINO_ACIDS);
		  }

	  }
	  else if( type.size() > 0)
	  {
		  std::cerr << "ERROR: no defined similarity type:<" << type << ">" << "\n";
		  exit( -1);
	  }
  }
  read.close();
  read.clear();
}



// is called by ReadScores( scores, score_and_profile_file) from sweet_align.cpp
// Score_and_Profile-file looks like: 10  ScaleSimilarity  KD.txt triangular 13 and can have several command-lines

void ReadScoresForList
(
		boost::shared_ptr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORES,
		const std::string &FILE,
		const std::set<char> &DEFINED_AMINO_ACIDS
)
{
  std::ifstream 
    read( FILE.c_str());

	std::string
		line;

	if( !read)
	{
		std::cerr <<"ERROR: The score_and_profile_file [" << FILE << "] was not found" << "\n";
		exit( -1);
	}


  std::string
    file,
    type,
    tmp,
    second_file;
  
  double
    weight;
  size_t
    count( 0);

	// read the score and profile file
  while( read)
  {
	  type.clear();
	  tmp.clear();

#ifdef SECURE
	  read >> tmp >> tmp;
	  if( tmp == "")
	  {
		  continue;
	  }
	  if( !mystr::IsNumerical( tmp) || ( weight = mystr::ConvertStringToNumericalValue< double>( tmp)) < 0)
	  {
		  std::cerr << "ERROR: weight has to be positive number or zero, found: <" << tmp << ">" << "\n";
		  exit( -1);
	  }
#else
	  read >> tmp >> weight;
#endif
	  read >> tmp >> type;


		// checks if the type is set to SequenceSimilartity
      if( type == "SequenceSimilarity")
      {

		  read >> tmp >> file;
		  DebugWrite( "SequenceSimilarity" << " " << file << " " << weight);

		  ScoreSequenceSimilarity score_seq_sim( file, DEFINED_AMINO_ACIDS);


		  if( SCORES)
		  {
			  boost::shared_ptr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
				  seq_sim( new ScoreSequenceSimilarity( score_seq_sim));
			  SCORES =  SCORES + weight * seq_sim;
		  }
		  else
		  {
			  SCORES = boost::shared_ptr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreSequenceSimilarity( score_seq_sim));
			  SCORES = weight * SCORES;
		  }
      }

		// checks for the other two possibilites: ProfileSimilarity and ScaleSimilarity
      else if( type == "ProfileSimilarity" || type == "ScaleSimilarity"|| type == "UniversalProfileSimilarity")
      {
    	  std::getline( read, line);
			// SCORES will be zero in the first loop
		  if( SCORES)
		  {
			  boost::shared_ptr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
			  seq_sim( new ScoreProfileSimilarity( count));
			  SCORES =  SCORES + weight * seq_sim;
		  }
		  else
		  {
			  SCORES = boost::shared_ptr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarity( count));
			  // SCORES is a pointer to a function object !

			  // that function object has as template specifications a pair of GeneralizedAA and a double
			  // the right side of the equation is calling the constructor of the shared pointer
			  // the constructor of the shared pointer that is called is the one that takes a 'new object' as input
			  // ptr< parent> name( new parent);
			  // now parent has a child class:
			  // ptr< parent> name( new child);
			  // child itself has a constructor child( value)
			  // ptr< parent> name( new child( value));

			  SCORES = weight * SCORES;

			  // now SCORES is still a pointer to a function object, but the function object is actually a sum function

		  }
		  ++count;
      }
      else if( type.size() > 0)
      {
    	  std::cout << "WARNING: no defined type found: " << type << "\n" << "\n";
      }
  }
  read.close();
  read.clear();
}



void ReadMultipleFileNames( std::vector< Sequence> &SEQUENCES, const std::string &FILE, std::set<char> &DEFINED_AMINO_ACIDS)
{
	  std::ifstream
	    read( FILE.c_str());

		if( !read)
		{
			std::cerr <<"ERROR: The file ["<< FILE << "] was not found" << "\n";
			exit( -1);
		}

	  std::string
	    line;

	  while( read)
	  {
		  line.clear();
		  read >> line;
		  if( line.size() > 0)
		  {
			  Sequence sequence( line);
			  sequence.SetFastaHeader( ReadSequence( sequence, line, DEFINED_AMINO_ACIDS));
			  SEQUENCES.push_back( sequence);
		  }
	  }
}


#endif /* SEQUENCE_READER_H */

