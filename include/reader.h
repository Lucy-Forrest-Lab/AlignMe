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


#ifndef READER_H
#define READER_H

#include <set>

#include "sequence.h"
#include "sum_function.t.h"

#include "../include/macro_functions_read_write.h"


//#include "sliding_window.h"
//#include "string_functions.h"
//
//
//#ifdef POSITION_SPECIFIC_SIMILARITY
//#include "score_position_specific_similarity.h"
//#else
//#include "sequence_functions.h"
//#endif



// the itr goes through the amino-acid-sequence that has been read in (f.e. first_sequence or second_sequence)
// itr->GetType () returns the m_Type (which is in this case an AA) that is stored at the position to which itr points to
// This value will be searched in Scale_Map in which values corresponding to a certain AA are stored as pairs
// with -> second you have access to this corresponding value
// This value will than be added to the sequence hat the position that itr points to
// Because of AddNewProfile(value) the value will be pushed back in m_Profiles of the corresponding GAA

void ReadScaleIntoSequence( Sequence &SEQUENCE,const std::map<char, double> &SCALE_MAP);


struct
CompareSecondOfPair
{
	bool operator()( const std::pair< char, double>& left, const std::pair< char, double>& right) const
	{
		return left.second < right.second;
	}
};


// called from  ReadSimilarityScoreFile
// with ReadScale( FIRST, file, sliding_window_type, window_size);
// with ReadScale( SECOND, file, sliding_window_type, window_size);
// file is the scale_file, f.e. KD.txt
void ReadScale
(
		Sequence &SEQUENCE,
		const std::string &SCALE_FILE,
		const std::string &SLIDING_WINDOW_TYPE,
		const int &WINDOW_SIZE,
		std::set<char> &DEFINED_AMINO_ACIDS,
		const float &THRESHOLD = 0.0
);



void ReadScaleIntoGappySequence
(
		Sequence &SEQUENCE,
		const std::string &SCALE_FILE,
		const std::string &SLIDING_WINDOW_TYPE,
		const size_t &WINDOW_SIZE,
		std::set<char> &DEFINED_AMINO_ACIDS
);


// Is called by ReadProfile( FIRST, first_file); or ReadProfile( SECOND, first_file)
// where FIRST and SECOND are constructors of GeneralizedAminoAcids
void  ReadProfile( Sequence &SEQUENCE, const std::string &PROFILE_FILE);



// Is called by ReadProfile( FIRST, first_file); or ReadProfile( SECOND, first_file)
// where FIRST and SECOND are constructors of GeneralizedAminoAcids
void  ReadUniversalProfile( Sequence &SEQUENCE, const std::string &PROFILE_FILE, const size_t &COLUMN, const int  &HEADER_LINES);




// Called from sweet_aling.cpp with   ReadSimilarityScoreFile( first_sequence, second_sequence, score_and_profile_file);
// Score_and_Profile-file looks like: 10  ScaleSimilarity  KD.txt triangular 13 and can have several command-lines
// Calculates a score for each position

void  ReadSimilarityScoreFile
(
		ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORES,
		Sequence &FIRST,
		Sequence &SECOND,
		const std::string &SCORE_AND_PROFILE_FILE,
		std::set<char> &DEFINED_AMINO_ACIDS,
		std::vector <std::string> &WRITE_PROFILES_HEADER1,
		std::vector <std::string> &WRITE_PROFILES_HEADER2,
		std::string &FIRST_FASTA_ID,
		std::string &SECOND_FASTA_ID
);




void
ReadSimilarityScoreFileMsa
(
		ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORES,
		std::vector< Sequence> &MSA1,
		std::vector< Sequence> &MSA2,
		const std::string &SCORE_AND_PROFILE_FILE,
		std::set<char> &DEFINED_AMINO_ACIDS
);




//////////////////////////////////////////////////////////////////////

// is called by ReadScores( scores, score_and_profile_file) from sweet_align.cpp
// Score_and_Profile-file looks like: 10  ScaleSimilarity  KD.txt triangular 13 and can have several command-lines





std::string
ReadSequence( Sequence &SEQUENCE, const std::string &FASTA_FILE, std::set<char> &DEFINED_AMINO_ACIDS);


void
ReadMultipleSequences
(
		std::vector< Sequence> &SEQUENCES,
		const std::string &FASTA_FILE,
		std::set<char> &DEFINED_AMINO_ACIDS
);


#endif /* READER_H */

