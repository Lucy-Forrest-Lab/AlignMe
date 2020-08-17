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


#include "../include/reader.h"
#include "../include/sliding_window.h"
#include "../include/string_functions.h"
#include "../include/score_sequence_similarity.h"
#include "../include/score_sequence_similarity_profile_dependent.h"
#include "../include/score_profile_similarity.h"
#include "../include/score_profile_similarity_linear_normalized.h"


#ifdef POSITION_SPECIFIC_SIMILARITY
#include "../include/score_position_specific_similarity.h"
#endif


// the itr goes through the amino-acid-sequence that has been read in (f.e. first_sequence or second_sequence)
// itr->GetType () returns the m_Type (which is in this case an AA) that is stored at the position to which itr points to
// This value will be searched in Scale_Map in which values corresponding to a certain AA are stored as pairs
// with -> second you have access to this corresponding value
// This value will than be added to the sequence hat the position that itr points to
// Because of AddNewProfile(value) the value will be pushed back in m_Profiles of the corresponding GAA

void ReadScaleIntoSequence(Sequence &SEQUENCE,const std::map<char, double> &SCALE_MAP)
{
	DebugWrite(__FUNCTION__);
	for(  std::vector< GeneralizedAminoAcid>::iterator itr = SEQUENCE.begin(); itr != SEQUENCE.end(); ++itr)
    {
		// std::map<>::find() returns an iterator to key/value std::pair // has to be used rather then operator [] because map is passed as const
		itr->AddNewProfile( SCALE_MAP.find( itr->GetType())->second);
		// (*itr).AddNewProfile( SCALE_MAP.find( (*itr).GetType()->second);
    }
}

//
//struct
//CompareSecondOfPair
//{
//	bool operator()( const std::pair< char, double>& left, const std::pair< char, double>& right) const
//	{
//		return left.second < right.second;
//	}
//};


void ReadBetaStrandScore
(
		Sequence &SEQUENCE,
		const std::string &BETASTRANDSCORE_FILE,
		const int &WINDOW_SIZE,
		std::set<char> &DEFINED_AMINO_ACIDS
)
{

	std::map<char, double>
		pore,
		membrane;
	double
		scale_value1,
		scale_value2;
	char
		aa_type;

	// read scale from file into map
	std::ifstream
		read(BETASTRANDSCORE_FILE.c_str());

	if( !read)
	{
		std::cerr << "ERROR: The scale file [" << BETASTRANDSCORE_FILE << "] was not found" << "\n";
		exit( -1);
	}

	while( read)
	{
		read >> aa_type >> scale_value1 >> scale_value2;
		pore.insert( std::make_pair( aa_type, scale_value1));
		membrane.insert( std::make_pair( aa_type, scale_value2));
	}

	ReadBetaStrandIntoSequence( SEQUENCE, pore, membrane, WINDOW_SIZE);
}


// called from ReadSimilarityScoreFile
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
		const float &THRESHOLD
)
{
	std::string
		aa_str,
		value_str;
	std::map<char, double>
		scale;
	double
		scale_value;
	char
		aa_type;

	// read scale from file into map
	std::ifstream
		read(SCALE_FILE.c_str());

	if( !read)
	{
		std::cerr << "ERROR: The scale file [" << SCALE_FILE << "] was not found" << "\n";
		exit( -1);
	}



	// read the scale_file and assign to every amino_acid (first column) the corresponding
	// hydrophobicity-value (second column)
	// store a pair (AA-Type, Hydrophocity) in scale

	while( read)
	{
//      read >> aa_type >> scale_value; // fast but unsafe !!

#ifdef SECURE

		aa_str.clear();
		read >> aa_str >> value_str;
		if( aa_str.size() == 0)
		{
			continue;
		}

		if( aa_str.size() > 1)
		{
			std::cerr << "ERROR: In scale files only single-letter amino acid keys are allowed (i.e. Y and not TYR)" << "\n";
			std::cerr << "ERROR: But found " << aa_str << " in scale [" << SCALE_FILE << "]"<< "\n";
			exit( -1);
		}

		aa_type = aa_str[0];

		if( !mystr::IsNumerical( value_str))
		{
			std::cerr << "ERROR: The entry " << aa_str << " in your scale file (" << SCALE_FILE << ") has a non-numeric value assigned, found: <" << value_str <<  ">!" << "\n";
			exit( -1);
		}

		scale_value = mystr::ConvertStringToNumericalValue< double>( value_str);

		if( int( aa_type) > 64 && int( aa_type) < 91)
		{
			if( scale.find( aa_type) != scale.end())
			{
				std::cout << "ERROR: The amino acid <" << aa_type << "> appears multiple times in the scale (" << SCALE_FILE << ") and therefore its definition is unclear. Please check for a correct formatting of the scale \n";
				exit(-1);
			}
			scale.insert( std::make_pair( aa_type, scale_value));
		}
		else
		{
			 //std::cout << "WARNING: the value "<< aa_type <<" from the scale " << SCALE_FILE << " does not appear in any of your sequences " << "\n";
		}
	}
	// check if the amino_acids that are in the sequences that have been read in, are in the scale
	// if an amino acid is in the sequence and not in the scale, than the program will quit
	for (std::set<char>::iterator acid_itr = DEFINED_AMINO_ACIDS.begin(); acid_itr!= DEFINED_AMINO_ACIDS.end(); ++acid_itr)
	{
		std::map<char, double>::const_iterator itr = scale.find(*acid_itr);
		if( itr == scale.end())
		{
			std::cerr << "ERROR: The amino acid " << *acid_itr<< " does not exist in the scale " << SCALE_FILE<< " that you want to use" << "\n";
			exit( -1);
		}
	}

#else
		read >> aa_type >> scale_value;
		scale.insert( std::make_pair( aa_type, scale_value));
	}
#endif



	// read scale into sequence with map, according to algorithm to use
	if( SLIDING_WINDOW_TYPE == "none")
    {
		ReadScaleIntoSequence( SEQUENCE, scale);
    }
	else if( SLIDING_WINDOW_TYPE == "rectangular")
	{
		RectangularSlidingWindow( SEQUENCE, scale, WINDOW_SIZE);
    }
	else if( SLIDING_WINDOW_TYPE == "triangular")
    {
		TriangularSlidingWindow( SEQUENCE, scale, WINDOW_SIZE);
    }
	else if( SLIDING_WINDOW_TYPE == "zigzag")
    {
		ZigZagSlidingWindow( SEQUENCE, scale, WINDOW_SIZE);
    }
	else if( SLIDING_WINDOW_TYPE == "sinusoidal")
    {
		SinoidSlidingWindow( SEQUENCE, scale, WINDOW_SIZE);
    }
//	else if( SLIDING_WINDOW_TYPE == "amphiphatic")
//	{
//		float delta = std::abs( std::max_element( scale.begin(), scale.end(), CompareSecondOfPair())->second - std::min_element( scale.begin(), scale.end(), CompareSecondOfPair())->second);
//		AmphiphaticSlidingWindow( SEQUENCE, scale, WINDOW_SIZE, delta, THRESHOLD);
//	}
//	else if( SLIDING_WINDOW_TYPE == "amphiphatic_derivatives")
//	{
//		float delta = std::abs( std::max_element( scale.begin(), scale.end(), CompareSecondOfPair())->second - std::min_element( scale.begin(), scale.end(), CompareSecondOfPair())->second);
//		AmphiphaticNextNeighborSlidingWindow( SEQUENCE, scale, WINDOW_SIZE, delta, THRESHOLD);
//	}
	else if( SLIDING_WINDOW_TYPE == "AmphiphilicityPure")
    {
		AmphiphilicityPure( SEQUENCE, scale, WINDOW_SIZE);
    }
	else if( SLIDING_WINDOW_TYPE == "AmphiphilicitySquared")
    {
		AmphiphilicitySquared( SEQUENCE, scale, WINDOW_SIZE);
    }
	else if( SLIDING_WINDOW_TYPE == "AmphiphilicityPureRectangular")
    {
		AmphiphilicityPureRectangular( SEQUENCE, scale, WINDOW_SIZE);
    }
	else if( SLIDING_WINDOW_TYPE == "AmphiphilicitySquaredRectangular")
    {
		AmphiphilicitySquaredRectangular( SEQUENCE, scale, WINDOW_SIZE);
    }
	else
	{
		std::cerr << "ERROR: Invalid window type: <" << SLIDING_WINDOW_TYPE << ">" << "\n";
		exit( -1);
	}

	read.close();
	read.clear();
}


void ReadScaleIntoGappySequence
(
		Sequence &SEQUENCE,
		const std::string &SCALE_FILE,
		const std::string &SLIDING_WINDOW_TYPE,
		const size_t &WINDOW_SIZE,
		std::set<char> &DEFINED_AMINO_ACIDS
)
{
	std::string
		aa_str,
		value_str;
	std::map<char, double> scale;
	double scale_value;
	char aa_type;

	// read scale from file into map
	std::ifstream read(SCALE_FILE.c_str());

	if( !read)
	{
		std::cerr << "ERROR: The scale file [" << SCALE_FILE << "] was not found" << "\n";
		exit( -1);
	}

	// read the scale_file and assign to every amino_acid (first column) the corresponding
	// hydrophobicity-value (second column)
	// store a pair (AA-Type, Hydrophobicity) in scale

	while( read)
    {

#ifdef SECURE

		aa_str.clear();

		read >> aa_str >> value_str;

		if( aa_str.size() == 0)
		{
			continue;
		}
		else if( aa_str.size() > 1)
		{
			std::cerr << "ERROR: Only single-letter amino acid keys with a corresponding value are allowed in scales but in the scale:  " << SCALE_FILE << "the value " << aa_str <<  "was found \n";
			exit( -1);
		}

		aa_type = aa_str[0];

		if( !mystr::IsNumerical( value_str))
		{
			std::cerr << "ERROR: The entry " << value_str << " in your provided scale file (" << SCALE_FILE << ") is not numeric!" << "\n";
			exit( -1);
		}

		scale_value = mystr::ConvertStringToNumericalValue< double>( value_str);

		if( ( int( aa_type) > 64 && int( aa_type) < 91) || aa_type == '-')
		{
			if( scale.find( aa_type) != scale.end())
			{
				std::cout << "ERROR: The amino acid <" << aa_type << "> appears multiple times in the scale (" << SCALE_FILE << ") and therefore its definition is unclear. Please check for a correct formatting of the scale \n";
				exit(-1);
			}
			scale.insert( std::make_pair( aa_type, scale_value));
		}
		else
		{
			//std::cout << "WARNING: the value "<< aa_type <<" from the scale " << SCALE_FILE << " does not appear in any of your sequences " << "\n";
		}
	}
	// check if the amino_acids that are in the sequences that have been read in, are in the scale
	// if an amino acid is in the sequence and not in the scale, than the program will quit
	for (std::set<char>::iterator acid_itr = DEFINED_AMINO_ACIDS.begin(); acid_itr!= DEFINED_AMINO_ACIDS.end(); ++acid_itr)
	{
		std::map<char, double>::const_iterator itr = scale.find(*acid_itr);
		if( *acid_itr != '-' && itr == scale.end())
		{
			std::cerr << "ERROR: The amino acid " << *acid_itr<< " does not exist in the scale " << SCALE_FILE<< " that you want to use" << "\n";
			exit( -1);
		}
	}

#else
		read >> aa_type >> scale_value;
		scale.insert( std::make_pair( aa_type, scale_value));
	}
#endif


	// read scale into sequence with map, according to algorithm to use
  if( SLIDING_WINDOW_TYPE == "none")
  {
      ReadScaleIntoSequence( SEQUENCE, scale);
  }
  else if (SLIDING_WINDOW_TYPE == "triangular_msa")
  {
	  TriangularSlidingWindowIgnoringGaps( SEQUENCE, scale, WINDOW_SIZE);
  }
  else
  {
	  std::cerr << "ERROR: For multiple sequence alignments the only window type that is currently allowed is 'triangular_msa'" << "\n";
	  exit( -1);
  }
  read.close();
  read.clear();
}

// Is called by ReadPositionSpecificSimilarityMatrix(SECOND, second_file, type) or 	 ReadPositionSpecificSimilarityMatrix(FIRST, first_file, type);
// where FIRST and SECOND are constructors of GeneralizedAminoAcids

void  ReadPositionSpecificSimilarityMatrix( Sequence &SEQUENCE, const std::string &PSSM_FILE, const std::string &TYPE)
{
	std::ifstream
		read( PSSM_FILE.c_str());

	std::string
		tmp,
		line;

// start is set to 2 because the first two entries in a row of a PSSM are <sequence position> <amino acid type>
// after this, there are the values for the amino acids
	int
		counter = 0,
		value,
		start = 2;

	size_t
		end;

	 std::vector< std::string>
		splitter;

	std::vector<char>
		pssm_acids;

// In contrast to the normal or profile PositionSpecificSimilarity, fdotf considers the raw counts and not the normalized ones.
	if (TYPE == "fdotfPositionSpecificSimilarity")
	{
		start = start+20;
	}
	end = start + 20;

	if( !read)
	{
		std::cerr << "ERROR: The PSSM file [" << PSSM_FILE << "] was not found" << "\n";
		exit( -1);
	}

	// skip 2 lines at the beginning of the file ( 1st: empty, 2nd: description)
	std::getline (read, line);
	std::getline (read, line);

	// 3rd line contains the amino acids written above the columns as headers
	std::getline (read, line);

	splitter = mystr::SplitString(line, " ");

#ifdef SECURE
	if (splitter.size() < 20 )
	{
		std::cerr << "ERROR: A list of all amino acids was not found the PSSM file [" << PSSM_FILE << "] \n";
		std::cerr << "This list should be in the third row of this file!\n";
		exit( -1);
	}
#endif

	for (size_t pssm_pos = 0; pssm_pos < 20; ++pssm_pos)
	{
		pssm_acids.push_back(splitter[pssm_pos][0]);
	}
	if (SEQUENCE.size() == 0)
	{

		while (read)
		{
			std::getline( read, line);
			if (line != "")
			{
				splitter = mystr::SplitString(line, " ");
			// skipping lines that are empty or contain calculated scores
				if ( splitter[0] != "" &&   splitter[0] != "K" &&   splitter[0] != "Standard" &&   splitter[0] != "PSI")
				{

					GeneralizedAminoAcid gaa;

					if (splitter.size() < end)
					{
						std::cerr << "The line < "<< line << " > in your PSSM is not formatted correctly. There have to be at least "<< end << " columns in a PSSM if alignments based on "<< TYPE << " should be generated\n";
						std::cerr << "Please check the PSSM file [" << PSSM_FILE << "] for a correct formatting \n";
						exit( -1);
					}

					for (size_t pssm_pos = start; pssm_pos < end; ++pssm_pos)
					{

#ifdef SECURE
						if( !mystr::IsNumerical(splitter[pssm_pos]))
						{
							std::cerr << "ERROR: The value: <" << value << "> is not numerical although it is expected to be \n";
							std::cerr << "Please check the PSSM file [" << PSSM_FILE << "] for a correct formatting \n";
							exit( -1);
						}
#endif
						value = mystr::ConvertStringToNumericalValue< double> (splitter[pssm_pos]);

						gaa.SetPSSM(pssm_acids[pssm_pos-2], value);
					}
					SEQUENCE.push_back( gaa);
				}
			}
		}
	}


	else
	{
		std::vector< GeneralizedAminoAcid>::iterator
			itr = SEQUENCE.begin();

		while( !read.eof() && itr != SEQUENCE.end())
		{
			if( itr->GetType() == '-')
			{
				continue;
			}

			std::getline( read, line);
			if (line != "")
			{
				splitter = mystr::SplitString(line, " ");
				// skipping lines that are empty or contain calculated scores

				if ( splitter[0] != "" &&   splitter[0] != "K" &&   splitter[0] != "Standard" &&   splitter[0] != "PSI")
				{

					if (splitter.size() < end)
					{
						std::cerr << "The line < "<< line << " > in your PSSM_file is not formatted correctly. There have to be at least "<< end << "columns in a PSSM if alignments based on ProfileSimilarity should be generated\n";
						std::cerr << "Please check the PSSM file [" << PSSM_FILE << "] for a correct formatting \n";
						exit( -1);
					}

					for (size_t pssm_pos = start; pssm_pos < end; ++pssm_pos)
					{
#ifdef SECURE
						if( !mystr::IsNumerical(splitter[pssm_pos]))
						{
							std::cerr << "ERROR: In the PSSM, in the line <"<< line << ">, the value: <" << splitter[pssm_pos] << "> is not numerical although it is expected to be \n";
							std::cerr << "Please check the PSSM file [" << PSSM_FILE << "] for a correct formatting! \n";
							exit( -1);
						}
#endif
						value = mystr::ConvertStringToNumericalValue< double> (splitter[pssm_pos]);

						itr->SetPSSM(pssm_acids[pssm_pos-2], value);
					}
					counter++;
					++itr;
				}
			}
		}

#ifdef SECURE

		if( itr != SEQUENCE.end())
		{
			std::cerr << "ERROR: The PSSM-profile obtained from the file < " << PSSM_FILE << " > is shorter than the corresponding sequence that you submitted! \n";
			std::cerr << "If you submit sequences and PSSMs at the same time, please be sure that they correspond to one another and that the formatting concurs with the guidelines mentioned in the manual \n";
			exit( -1);
		}

		while( !read.eof())
		{
			std::getline( read, line);
			if (line != "")
			{
				splitter = mystr::SplitString(line, " ");
				// skipping lines that are empty or contain calculated scores

				if ( splitter[0] != "" &&   splitter[0] != "K" &&   splitter[0] != "Standard" &&   splitter[0] != "PSI")
				{
					std::cerr << "ERROR: The PSSM-profile obtained from the file < " << PSSM_FILE << " > is longer than the corresponding sequence that you submitted! \n";
					std::cerr << "If you submit sequences and PSSMs at the same time, please be sure that they correspond to one another and that the formatting concurs with the guidelines mentioned in the manual \n";
					exit( -1);
				}
			}
		}

#endif
	}
	read.close();
	read.clear();

}



// Is called by ReadUniversalProfile( FIRST, first_file); or ReadUniversalProfile( SECOND, first_file)
// where FIRST and SECOND are constructors of GeneralizedAminoAcids
void  ReadUniversalProfile( Sequence &SEQUENCE, const std::string &PROFILE_FILE, const size_t &COLUMN, const int  &HEADER_LINES)
{

	std::ifstream
		read( PROFILE_FILE.c_str());

	std::string
		tmp,
		line;

	double
		value;

	 std::vector< std::string>
		splitter;


	if( !read)
	{
		if ( PROFILE_FILE != "")
		{
			std::cerr << "ERROR: The profile file [" << PROFILE_FILE << "] was not found" << "\n";
		}
		else
		{
			std::cerr << "ERROR: You have to provide 2 profile files to make an alignment based on similarity of profiles" << "\n";
		}
		exit( -1);
	}

	for (int i = 0; i != HEADER_LINES; i++)
	{
		std::getline( read, line);
	}
	// If the GeneralizedAminoAcidSequence is empty an nothing has been added to it yet
	// then the profile_file will be opened and as long as there are values been read in:
	// m_Type (X) will be added because of GeneralizedAminoAcid()
	// m_Profiles (value) will also be added at this position

	if (SEQUENCE.size() == 0)
	{
		while (read)
		{
			value = std::numeric_limits<double>::max();
			std::getline( read, line);

			if( line != "")
			{
				splitter = mystr::SplitString(line, " ");

#ifdef SECURE
				if (splitter.size() < COLUMN)
				{
					std::cerr << "ERROR: Your profile in profile file (" << PROFILE_FILE << ") does not contain enough columns \n";
					std::cerr << "COLUMN "<< COLUMN <<" can not be accessed because only "<< splitter.size() <<" columns are available  \n";
					exit( -1);
				}
				if( !mystr::IsNumerical(splitter[COLUMN-1]))
				{
					std::cerr << "ERROR: The value: <" << splitter[COLUMN-1] << "> in column " << COLUMN << " of profile file (" << PROFILE_FILE << ") is not numerical" << "\n";
					exit( -1);
				}
#endif
				value = mystr::ConvertStringToNumericalValue< double>( splitter[COLUMN-1]);
				SEQUENCE.push_back(GeneralizedAminoAcid());
				SEQUENCE.back().AddNewProfile(value);
			}
		}
	}

	// If there are already values stored in the GAA-Sequence:
	// the itr goes through the whole sequence and points to a certain GAA
	// a value corresponding to the profile_file will be added to this certain GAA in m_Profiles
	// m_Profiles can contain a list of different values corresponding to a certain GAA
	else
	{
//		std::cout <<  __FUNCTION__ << " sequence given" << "\n";
		std::vector< GeneralizedAminoAcid>::iterator
			itr = SEQUENCE.begin();
		while( !read.eof() && itr != SEQUENCE.end())
		{
			if( itr->GetType() == '-')
			{
				continue;
			}

			value = std::numeric_limits<double>::max();
			std::getline( read, line);

			if( line != "")
			{
				splitter = mystr::SplitString(line, " ");

#ifdef SECURE
				if (splitter.size() < COLUMN)
				{
					std::cerr << "ERROR: Your profile in profile file (" << PROFILE_FILE << ") does not contain enough columns \n";
					std::cerr << "COLUMN "<< COLUMN <<" can not be accessed because only "<< splitter.size() <<" columns are available  \n";
					exit( -1);
				}
				if( !mystr::IsNumerical(splitter[COLUMN-1]))
				{
					std::cerr << "ERROR: The value: <" << splitter[COLUMN-1] << "> in column " << COLUMN << " of profile file (" << PROFILE_FILE << ") is not numerical" << "\n";
					exit( -1);
				}
#endif
				value = mystr::ConvertStringToNumericalValue< double>( splitter[COLUMN-1]);
				itr->AddNewProfile( value);
				++itr;
			}
		}

		read >> tmp;

#ifdef SECURE
		if( !read.eof())
		{
			std::cerr << "ERROR: The profile in " << PROFILE_FILE << " is longer than the corresponding sequence that you submitted!" << "\n";
			std::cerr << "If you submit sequences and profiles at the same time, please be sure that they correspond to one another" << "\n";
			exit( -1);
		}
		if( itr != SEQUENCE.end())
		{
			std::cerr << "ERROR: The profile in " << PROFILE_FILE << " is shorter than the corresponding sequence that you submitted!" << "\n";
			std::cerr << "If you submit sequences and profiles at the same time, please be sure that they correspond to one another" << "\n";
			exit( -1);
		}
#endif
	}
	read.close();
	read.clear();
}




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
)
{
	DebugWrite(  __FUNCTION__ );
 std::string
	 line,
	 file,
	 type,
	 first_file,
	 second_file,
	 tmp,
	 sliding_window_type;


	 std::vector< std::string>
		splitter;

  double
	  threshold = 0.5,   /// for amphiphatic window algorithms
	  weight;

  int
	  splitter_size,
	  window_size,
	  header_lines;


  size_t
	column,
    count( 0);


  bool
	  file_not_empty = "FALSE";

  std::ifstream 
    test_read,
    read( SCORE_AND_PROFILE_FILE.c_str());
  

	if( !read)
	{
		std::cerr << "ERROR: The score_and_profile_file [" << SCORE_AND_PROFILE_FILE << "] was not found" << "\n";
		exit( -1);
	}

  while (read)
  {
	  std::getline( read, line);

	  if (line != "")
	  {
		  file_not_empty = "TRUE";
		  splitter = mystr::SplitString(line, " ");
		  splitter_size = splitter.size();

#ifdef SECURE
		  if (splitter_size < 6)
		  {
			  std::cerr << "ERROR: A line in your similarity_score_file contains too less parameters! Please check if you have use the correct file and have a look at the manual on how to set up this file correctly!\n";
			  exit( -1);
		  }
#endif
		  weight = mystr::ConvertStringToNumericalValue< double>(splitter[1]);

#ifdef SECURE
		  if( !mystr::IsNumerical( splitter[1]) || (weight  < 0))
		  {
			  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the weight has to be a positive number, or zero, found: <" << splitter[1] << ">" << "\n";
			  exit( -1);
		  }
#endif
		  type = splitter[3];
		  if( type == "SequenceSimilarity" )
		  {
#ifdef SECURE
			  if (splitter_size != 6)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" file: <filename>\n";
				  exit( -1);
			  }

			  test_read.open( splitter[5].c_str());

			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the substitution matrix [" << splitter[5] << "] was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();
#endif
			  std::map< std::string, double>
				  map = ReadSubstitutionMatrix( splitter[5] , DEFINED_AMINO_ACIDS);

			  if( SCORES)
			  {
				  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
					  seq_sim( new ScoreSequenceSimilarity( map));
				  SCORES =  SCORES + weight * seq_sim;
			  }
			  else
			  {
				  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreSequenceSimilarity( map));
				  SCORES = weight * SCORES;
			  }


		  }
		  else if ( type == "LinearNormalizedSequenceSimilarity")
		  {
#ifdef SECURE
			  if (splitter_size != 8)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" file: <filename> min_cutoff <double> \n";
				  exit( -1);
			  }

			  test_read.open( splitter[5].c_str());

			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the substitution matrix [" << splitter[5] << "] was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();

			  if( !mystr::IsNumerical( splitter[7]))
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the value for min_cutoff has to be a number, found <" << splitter[7] << ">" << "\n";
				  exit( -1);
			  }
#endif
			  std::map< std::string, double>
				  map = ReadSubstitutionMatrix( splitter[5] , DEFINED_AMINO_ACIDS);

			 double min_cutoff = mystr::ConvertStringToNumericalValue<int>( splitter[7]);

			  map = ScaleSubstitutionMatrixToUnity( map, min_cutoff );

			  if( SCORES)
			  {
				  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
					  seq_sim( new ScoreSequenceSimilarity( map));
				  SCORES =  SCORES + weight * seq_sim;
			  }
			  else
			  {
				  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreSequenceSimilarity( map));
				  SCORES = weight * SCORES;
			  }



		  }
		  else if ( type == "NormalizedProfileDependentSubstitutionMatrix")
		  {
#ifdef SECURE
			  if (splitter_size != 8)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" matrix1: <filename> matrix2: <filename> profile: <int> threshold: <double> \n";
				  exit( -1);
			  }

			  test_read.open( splitter[5].c_str());

			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the substitution matrix [" << splitter[5] << "] was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();
			  test_read.open( splitter[7].c_str());

			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the substitution matrix [" << splitter[7] << "] was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();
			  // profile_id
			  if( !mystr::IsNumerical( splitter[9]))
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the profile which should be used for the assignment has to be referred to as a number, found <" << splitter[9] << ">" << "\n";
				  exit( -1);
			  }

#endif
			 double profile_id = mystr::ConvertStringToNumericalValue<int>( splitter[9]);
			  // threshold
#ifdef SECURE
			  if( !mystr::IsNumerical( splitter[11]))
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the threshold that should be used for the assignment of different matrices on a profile has to be a number,  found <" << splitter[11] << ">" << "\n";
				  exit( -1);
			  }
#endif
			  double threshold_on_profile = mystr::ConvertStringToNumericalValue<double>( splitter[11]);
			  std::map< std::string, double>
				  upper,
				  lower;

			  upper = ReadSubstitutionMatrix( splitter[5], DEFINED_AMINO_ACIDS);
			  lower = ReadSubstitutionMatrix( splitter[7], DEFINED_AMINO_ACIDS);

			  if( SCORES)
			  {
				  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
					  seq_sim( new ScoreProfileDependentSequenceSimilarity( upper, lower,  profile_id,  threshold_on_profile));
				  SCORES =  SCORES + weight * seq_sim;
			  }
			  else
			  {
				  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileDependentSequenceSimilarity( upper, lower, profile_id,  threshold_on_profile));
				  SCORES = weight * SCORES;
			  }


		  }

		  else if( type == "PositionSpecificSimilarity" || type == "MinPositionSpecificSimilarity"  || type == "ProfilePositionSpecificSimilarity" || type == "fdotfPositionSpecificSimilarity")
		  {
#ifdef SECURE
			  if (splitter_size != 8)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" file: <filename> min_cutoff <double> \n";
				  exit( -1);
			  }
#endif
			  first_file = splitter[5];
			  second_file = splitter[7];
			  test_read.open(  first_file.c_str());

#ifdef SECURE
			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the substitution matrix [" <<   first_file << "] was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();
			  test_read.open( second_file.c_str());

			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the substitution matrix [" << second_file << "] was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();
#endif
			  ReadPositionSpecificSimilarityMatrix(FIRST, first_file, type);
			  ReadPositionSpecificSimilarityMatrix(SECOND, second_file, type);


			  if( type == "PositionSpecificSimilarity")
			  {
				  if( SCORES)
				  {
					  SCORES =  SCORES + weight * ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScorePositionSpecificSimilarity());
				  }
				  else
				  {
					  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScorePositionSpecificSimilarity());
					  SCORES = weight * SCORES;
				  }
			  }

			  else if( type == "MinPositionSpecificSimilarity")
		      {

				  if( SCORES)
				  {
					  SCORES =  SCORES + weight * ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new MinScorePositionSpecificSimilarity());
				  }
				  else
				  {
					  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new MinScorePositionSpecificSimilarity());
					  SCORES = weight * SCORES;
				  }
		      }

		      else if( type == "fdotfPositionSpecificSimilarity")
		      {

				  if( SCORES)
				  {
					  SCORES =  SCORES + weight * ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new fdotfScorePositionSpecificSimilarity());
				  }
				  else
				  {
					  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new fdotfScorePositionSpecificSimilarity());
					  SCORES = weight * SCORES;
				  }
		      }

		      else if( type == "ProfilePositionSpecificSimilarity")
		      {

				  if( SCORES)
				  {
					  SCORES =  SCORES + weight * ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ProfileScorePositionSpecificSimilarity());
				  }
				  else
				  {
					  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ProfileScorePositionSpecificSimilarity());
					  SCORES = weight * SCORES;
				  }
		      }
		  }

		  else if( type == "UniversalProfileSimilarity" || type == "LinearNormalizedProfileSimilarity" )
		  {
#ifdef SECURE
			  if (splitter_size != 12 &&  type == "UniversalProfileSimilarity")
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" file: <filename> min_cutoff <double> \n";
				  exit( -1);
			  }
			  else  if (splitter_size != 14 &&  type == "LinearNormalizedProfileSimilarity")
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" file: <filename> min_cutoff <double> \n";
				  exit( -1);
			  }
#endif
			  column = mystr::ConvertStringToNumericalValue< size_t>(splitter[5]);
#ifdef SECURE
			  if( !mystr::IsNumerical( splitter[5]) || ( column  < 0) )
			  {
				  std::cerr<< "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the column value for UniversalProfileSimilarity has to be a positive number or zero, found: <"<< tmp<< "> in " <<  SCORE_AND_PROFILE_FILE << "\n";
				  exit( -1);
			  }
#endif
			  header_lines = mystr::ConvertStringToNumericalValue< size_t>( splitter[7]);
#ifdef SECURE
			  if( !mystr::IsNumerical( splitter[7]) || ( header_lines < 0) )
			  {
				  std::cerr<< "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the header_lines value for UniversalProfileSimilarity has to be a positive number or zero, found: <"<< tmp<< "> in " <<  SCORE_AND_PROFILE_FILE << "\n";
				  exit( -1);
			  }
#endif

			  first_file = splitter[9];
			  second_file = splitter[11];
#ifdef SECURE
			  test_read.open(  first_file.c_str());

			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the first profile [" <<   first_file << "] was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();
			  test_read.open( second_file.c_str());

			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the second profile [" << second_file << "] was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();
#endif
			  ReadUniversalProfile( FIRST, first_file, column, header_lines);
			  ReadUniversalProfile( SECOND, second_file, column, header_lines);

			  WRITE_PROFILES_HEADER1.push_back(" based on profile " + first_file + ", extracted column: "+mystr::NumericalValueToString(column)+"\n");
			  WRITE_PROFILES_HEADER2.push_back(" based on profile " + second_file + ", extracted column: "+mystr::NumericalValueToString(column)+"\n");

			  if( type == "UniversalProfileSimilarity")
			  {
					// SCORES will be zero in the first loop
				  if( SCORES)
				  {
					  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
						  seq_sim( new ScoreProfileSimilarity( count));
					  SCORES =  SCORES + weight * seq_sim;
				  }
				  else
				  {
					  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarity( count));
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
			  else if (type == "LinearNormalizedProfileSimilarity" )
			  {
				  double
					  max_delta;
				  if( !mystr::IsNumerical(splitter[13]) || (max_delta = mystr::ConvertStringToNumericalValue<int>(splitter[13]) < 0))
				  {
					  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the value for max_delta has to be a positive number, or zero, found: <" << splitter[13] << ">" << "\n";
					  exit( -1);
				  }

				  if( SCORES)
				  {
					  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
						  seq_sim( new ScoreProfileSimilarityLinearNormalized( count, max_delta));
					  SCORES =  SCORES + weight * seq_sim;
				  }
				  else
				  {
					  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarityLinearNormalized( count, max_delta));

					  SCORES = weight * SCORES;
				  }
				  ++count;
			  }

		  }



		  else if( type == "BetaStrandScore" || type == "LinearNormalizedBetaStrandScore")
		  {
#ifdef SECURE
			  if (splitter_size != 8)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" file: <filename> min_cutoff <double> \n";
				  exit( -1);
			  }
#endif
			  file = splitter[5];
#ifdef SECURE
			  test_read.open(  file.c_str());

			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the file for the betastrand score  [" <<   file << "] was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();
#endif
			 tmp = splitter[7];
	  		 if( !mystr::IsNumerical( tmp) || ( window_size = mystr::ConvertStringToNumericalValue<int>(tmp)) < 0)
		  	 {
		  		  std::cerr<< "ERROR concerning the similarity_score_file: windowsize value for the BetaStrandScore has to be a positive number or zero, found: <"<< tmp << ">" << "\n";
		  		  exit( -1);
		  	 }
#ifdef SECURE
		  	 if( int(FIRST.size()) <= window_size || int(SECOND.size()) <= window_size)
		  	 {
		  		  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], sequences have to be at least as long as the size of the window requested! "<< "\n";			  std::cerr << "ERROR: Length of provided sequences, first: <" << FIRST.size() << ">, second: <" << SECOND.size() << ">" << "\n";
		  		  std::cerr << "Length of provided sequences, first: <" << FIRST.size() << ">, second: <" << SECOND.size() << ">" << "\n";
		  		  exit( -1);
		  	 }
#endif
		  	 ReadBetaStrandScore( FIRST, file, window_size, DEFINED_AMINO_ACIDS);
		  	 ReadBetaStrandScore( SECOND, file, window_size, DEFINED_AMINO_ACIDS);

	  		  WRITE_PROFILES_HEADER1.push_back(" based on beta-strand-values of " +  file + " for sequence "+  FIRST_FASTA_ID +  " with size " + mystr::NumericalValueToString(window_size)+"\n");
	  		  WRITE_PROFILES_HEADER2.push_back(" based on beta-strand-values of " +  file + " for sequence "+  SECOND_FASTA_ID + " with size " + mystr::NumericalValueToString(window_size)+"\n");


	  		if( type == "BetaStrandScore")
  			{
	  		  if( SCORES)
	  		  {
	  			  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
	  				  seq_sim( new ScoreProfileSimilarity( count));
	  			  SCORES =  SCORES + weight * seq_sim;
	  		  }
	  		  else
	  		  {
	  			  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarity( count));
	  			  SCORES = weight * SCORES;
	  		  }
	  		  ++count;
  			}
	  		else if (type == "LinearNormalizedBetaStrandScore")
  			{
  			  double
					  value1, value2, max_delta,
					  min( std::numeric_limits< double>::max()),
					  max( -std::numeric_limits< double>::max());
  			  tmp.clear();


  			  std::ifstream in_scale( file.c_str());

  			  while( in_scale)
  			  {
  				  tmp.clear();
  				  in_scale >> tmp >> value1 >> value2;
  				  if( tmp.empty())
  				  {
  					  continue;
  				  }
  				  min = std::min( value1, min);
  				  min = std::min( value2, min);
  				  max = std::max( value1, max);
  				  max = std::max( value2, max);
  			  }
  			  max_delta = fabs( max - min);
  			  in_scale.close();
  			  in_scale.clear();

  			  if( SCORES)
  			  {
  				  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
  					  seq_sim( new ScoreProfileSimilarityLinearNormalized( count, max_delta));
  				  SCORES =  SCORES + weight * seq_sim;
  			  }
  			  else
  			  {
  				  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarityLinearNormalized( count, max_delta));

  				  SCORES = weight * SCORES;
  			  }
  			  ++count;
  			}

		  }

		  else if( type == "ScaleSimilarity" || type == "LinearNormalizedScaleSimilarity")
		  {
#ifdef SECURE
			  if (splitter_size != 10)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" file: <filename> min_cutoff <double> \n";
				  exit( -1);
			  }
#endif
			  file = splitter[5];
#ifdef SECURE
			  test_read.open(  file.c_str());

			  if( !test_read)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the file called   [" <<   file << "] for the type "<< type << " was not found";
				  exit( -1);
			  }
			  test_read.close();
			  test_read.clear();
#endif
			  sliding_window_type = splitter[7];

#ifdef SECURE
			  if( FIRST.size() == 0 || SECOND.size() == 0)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], you have to provide two sequences if you want to make an alignment based on scales"<<"\n";
				  std::cerr << "Length of provided sequences, first sequence: <" << FIRST.size() << ">, second sequence: <" << SECOND.size() << ">" << "\n";
				  exit( -1);
			  }
#endif
			  window_size = mystr::ConvertStringToNumericalValue<int>(splitter[9]);

#ifdef SECURE
			  if( !mystr::IsNumerical(splitter[9]) || ( window_size < 0) )
			  {
				  std::cerr<< "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], the windowsize value for UniveralProfileSimilarity has to be a positive number or zero, found: <"<< tmp << ">" << "\n";
				  exit( -1);
			  }
			  if( int(FIRST.size()) <= window_size || int(SECOND.size()) <= window_size)
			  {
				  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], sequences have to be at least as long as the size of the window requested! "<< "\n";			  std::cerr << "ERROR: Length of provided sequences, first: <" << FIRST.size() << ">, second: <" << SECOND.size() << ">" << "\n";
				  std::cerr << " The length of provided sequences, first: <" << FIRST.size() << ">, second: <" << SECOND.size() << ">" << "\n";
				  exit( -1);
			  }
#endif

//			  if( sliding_window_type.find("amphiphatic") != std::string::npos){
//				  read >> tmp >> threshold;
//			  }

			  ReadScale( FIRST, file, sliding_window_type, window_size, DEFINED_AMINO_ACIDS, threshold);
			  ReadScale( SECOND, file, sliding_window_type, window_size, DEFINED_AMINO_ACIDS, threshold);

			  WRITE_PROFILES_HEADER1.push_back(" based on scale " +  file + " for sequence "+  FIRST_FASTA_ID + " using windowtype " + sliding_window_type + " with size " + mystr::NumericalValueToString(window_size)+"\n");
			  WRITE_PROFILES_HEADER2.push_back(" based on scale " +  file + " for sequence "+  SECOND_FASTA_ID + " using windowtype "+ sliding_window_type + " with size " + mystr::NumericalValueToString(window_size)+"\n");

			  if( type == "ScaleSimilarity")
			  {
				  if( SCORES)
				  {
					  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
						  seq_sim( new ScoreProfileSimilarity( count));
					  SCORES =  SCORES + weight * seq_sim;
				  }
				  else
				  {
					  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarity( count));
					  SCORES = weight * SCORES;
				  }
				  ++count;
			  }
			  else if( type == "LinearNormalizedScaleSimilarity" )
			  {

    			  double
					  max_delta,
					  value,
					  min( std::numeric_limits< double>::max()),
					  max( -std::numeric_limits< double>::max());
    			  tmp.clear();


    			  std::ifstream in_scale( file.c_str());

    			  while( in_scale)
    			  {
    				  tmp.clear();
    				  in_scale >> tmp >> value;
    				  if( tmp.empty())
    				  {
    					  continue;
    				  }
    				  min = std::min( value, min);
    				  max = std::max( value, max);
    			  }
    			  max_delta = fabs( max - min);
    			  in_scale.close();
    			  in_scale.clear();

    			  if( SCORES)
    			  {
    				  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
    					  seq_sim( new ScoreProfileSimilarityLinearNormalized( count, max_delta));
    				  SCORES =  SCORES + weight * seq_sim;
    			  }
    			  else
    			  {
    				  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarityLinearNormalized( count, max_delta));

    				  SCORES = weight * SCORES;
    			  }
    			  ++count;
			  }

		  }
		  else
		  {
			  std::cerr << "ERROR: In the similarity_score_file  [" << SCORE_AND_PROFILE_FILE << "], you provided the type "<< type <<" which is not supported by AlignMe. Please use a supported type! \n";
			  exit( -1);
		  }

	 }

  }
	read.close();
	read.clear();

	if (!file_not_empty)
	{
		  std::cerr << "ERROR: Your similarity_score_file: <" << SCORE_AND_PROFILE_FILE << "> is empty or contains an empty line at the beginning" << "\n";
		  exit( -1);
	}
}


/* now included in ReadSimilarityScoreFile
void ReadScores
(
		ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORES,
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
	  if(tmp == "none" || tmp == "")
	  {
		  continue;
	  }
	  if( tmp != "" && ( !mystr::IsNumerical( tmp) || ( weight = mystr::ConvertStringToNumericalValue< double>( tmp)) < 0))
	  {
		  std::cerr << "ERROR: weight has to be positive number or zero, found: <" << tmp << ">\n";
		  std::cerr << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		  exit( -1);
	  }
#else
	  read >> tmp >> weight;
#endif
	  read >> tmp >> type;


		// checks if the type is set to SequenceSimilartity
      if( type == "SequenceSimilarity" || type == "LinearNormalizedSequenceSimilarity")
      {
    	  std::getline( read, line);

       	  std::vector< std::string>
			  bits = mystr::SplitString( line);


		  file =  bits[1];
		  DebugWrite( type << " " << file << " " << weight);

		  std::map< std::string, double>
			  map = ReadSubstitutionMatrix( file, DEFINED_AMINO_ACIDS);

		  if( type == "LinearNormalizedSequenceSimilarity")
		  {
			  double
				  min_cutoff = std::numeric_limits< double>::max();
			  if (bits.size() == 4)
			  {
		       	  std::vector< std::string>
					  bits = mystr::SplitString( line);

		       	  tmp = bits[3];

				  if( tmp == "min_cutoff:")
				  {
					  min_cutoff = mystr::ConvertStringToNumericalValue< size_t>( bits[4]);
				  }
			  }

			  map = ScaleSubstitutionMatrixToUnity( map, min_cutoff);
		  }

		  if( SCORES)
		  {
			  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
				  seq_sim( new ScoreSequenceSimilarity( map));
			  SCORES =  SCORES + weight * seq_sim;
		  }
		  else
		  {
			  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreSequenceSimilarity( map));
			  SCORES = weight * SCORES;
		  }
      }
      else if( type == "NormalizedProfileDependentSubstitutionMatrix")
      {
    	  int profile_id;
    	  double threshold;
		  read >> tmp >> file;
		  std::map< std::string, double>
			  upper = ReadSubstitutionMatrix( file, DEFINED_AMINO_ACIDS);
		  read >> tmp >> file;
		  std::map< std::string, double>
			  lower = ReadSubstitutionMatrix( file, DEFINED_AMINO_ACIDS);

		  upper = ScaleSubstitutionMatrixToUnity( upper);
    	  lower = ScaleSubstitutionMatrixToUnity( lower);

		  read >> tmp >> profile_id >> tmp >> threshold;


		  if( SCORES)
		  {
			  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
				  seq_sim( new ScoreProfileDependentSequenceSimilarity( upper, lower, profile_id, threshold));
			  SCORES =  SCORES + weight * seq_sim;
		  }
		  else
		  {
			  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileDependentSequenceSimilarity( upper, lower, profile_id, threshold));
			  SCORES = weight * SCORES;
		  }

      }
#ifdef POSITION_SPECIFIC_SIMILARITY
      else if( type == "PositionSpecificSimilarity")
      {

    	  read >> tmp >> file >> tmp >> second_file;
//		  DebugWrite( "PositionSpecificSimilarity" << " " << file << " " << second_file << " " << weight);
//
//		  if( SCORES)
//		  {
//			  SCORES =  SCORES + weight * ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScorePositionSpecificSimilarity());
//		  }
//		  else
//		  {
//			  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScorePositionSpecificSimilarity());
//			  SCORES = weight * SCORES;
//		  }
      }

      else if( type == "MinPositionSpecificSimilarity")
      {

		  read >> tmp >> file >> tmp >> second_file;
		  DebugWrite( "MinPositionSpecificSimilarity" << " " << file << " " << second_file << " " << weight);

		  if( SCORES)
		  {
			  SCORES =  SCORES + weight * ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new MinScorePositionSpecificSimilarity());
		  }
		  else
		  {
			  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new MinScorePositionSpecificSimilarity());
			  SCORES = weight * SCORES;
		  }
      }

      else if( type == "fdotfPositionSpecificSimilarity")
      {

		  read >> tmp >> file >> tmp >> second_file;
		  DebugWrite( "MinPositionSpecificSimilarity" << " " << file << " " << second_file << " " << weight);

		  if( SCORES)
		  {
			  SCORES =  SCORES + weight * ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new fdotfScorePositionSpecificSimilarity());
		  }
		  else
		  {
			  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new fdotfScorePositionSpecificSimilarity());
			  SCORES = weight * SCORES;
		  }
      }

      else if( type == "ProfilePositionSpecificSimilarity")
      {

		  read >> tmp >> file >> tmp >> second_file;
		  DebugWrite( "ProfilePositionSpecificSimilarity" << " " << file << " " << second_file << " " << weight);

		  if( SCORES)
		  {
			  SCORES =  SCORES + weight * ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ProfileScorePositionSpecificSimilarity());
		  }
		  else
		  {
			  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ProfileScorePositionSpecificSimilarity());
			  SCORES = weight * SCORES;
		  }
      }
#endif
      // checks for the other two possibilites: ProfileSimilarity and ScaleSimilarity
      else if( type == "ProfileSimilarity" || type == "ScaleSimilarity"|| type == "UniversalProfileSimilarity" || type == "BetaStrandScore")
      {
    	  std::getline( read, line);
//			// SCORES will be zero in the first loop
//		  if( SCORES)
//		  {
//			  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
//				  seq_sim( new ScoreProfileSimilarity( count));
//			  SCORES =  SCORES + weight * seq_sim;
//		  }
//		  else
//		  {
//			  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarity( count));
//			  // SCORES is a pointer to a function object !
//
//			  // that function object has as template specifications a pair of GeneralizedAA and a double
//			  // the right side of the equation is calling the constructor of the shared pointer
//			  // the constructor of the shared pointer that is called is the one that takes a 'new object' as input
//			  // ptr< parent> name( new parent);
//			  // now parent has a child class:
//			  // ptr< parent> name( new child);
//			  // child itself has a constructor child( value)
//			  // ptr< parent> name( new child( value));
//
//			  SCORES = weight * SCORES;
//
//			  // now SCORES is still a pointer to a function object, but the function object is actually a sum function
//
//		  }
//		  ++count;
      }
      else if( type == "LinearNormalizedScaleSimilarity" || type == "LinearNormalizedProfileSimilarity" ||  type == "LinearNormalizedBetaStrandScore")
      {
    	  double max_delta;
       	  std::getline( read, line);
       	  std::vector< std::string>
			  bits = mystr::SplitString( line);

#ifdef SECURE
       	  if (( type == "LinearNormalizedScaleSimilarity" && bits.size() == 8) || ( type == "LinearNormalizedProfileSimilarity" && bits.size() == 10) )
       	  {
           	  if( bits[ bits.size() - 2] != "max_delta:")
           	  {
				  std::cerr << "ERROR: If you want to define a max_delta value than the correct flag is'max_delta:'\n";
				  exit( -1);
           	  }
			  double tmp_nr;
			  if( !mystr::IsNumerical(bits[ bits.size() - 1]) || (tmp_nr = mystr::ConvertStringToNumericalValue<int>(bits[ bits.size() - 1]) < 0))
			  {
				  std::cerr << "ERROR: The value for max_delta has to be a positive number, or zero, found: <" << bits[ bits.size() - 1] << ">" << "\n";
				  exit( -1);
			  }

       	  }
#endif
       	  if( bits[ bits.size() - 2] == "max_delta:")
    	  {
    		  max_delta = mystr::ConvertStringToNumericalValue< double>( bits[ bits.size() - 1]);
    	  }
    	  else
    	  {
    		  // determine max delta from scale
    		  if( type == "LinearNormalizedScaleSimilarity")
    		  {
    			  double
					  value,
					  min( std::numeric_limits< double>::max()),
					  max( -std::numeric_limits< double>::max());
    			  tmp.clear();


    			  std::ifstream in_scale( bits[1].c_str());

    			  while( in_scale)
    			  {
    				  tmp.clear();
    				  in_scale >> tmp >> value;
    				  if( tmp.empty())
    				  {
    					  continue;
    				  }
    				  min = std::min( value, min);
    				  max = std::max( value, max);
    			  }
    			  max_delta = fabs( max - min);
    			  in_scale.close();
    			  in_scale.clear();
    		  }
    		  else if ( type == "LinearNormalizedBetaStrandScore")
    		  {
    			  double
					  value1, value2,
					  min( std::numeric_limits< double>::max()),
					  max( -std::numeric_limits< double>::max());
    			  tmp.clear();


    			  std::ifstream in_scale( bits[1].c_str());

    			  while( in_scale)
    			  {
    				  tmp.clear();
    				  in_scale >> tmp >> value1 >> value2;
    				  if( tmp.empty())
    				  {
    					  continue;
    				  }
    				  min = std::min( value1, min);
    				  min = std::min( value2, min);
    				  max = std::max( value1, max);
    				  max = std::max( value2, max);
    			  }
    			  max_delta = fabs( max - min);
    			  in_scale.close();
    			  in_scale.clear();
    		  }
#ifdef SECURE
    		  else
    		  {
    			  std::cout << "ERROR: for linear normalized profiles a 'max_delta:' value has to be provided in the similarity score file" << std::endl;
    			  exit( -1);
    		  }
#endif

    	  }
 			// SCORES will be zero in the first loop
		  if( SCORES)
		  {
			  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
				  seq_sim( new ScoreProfileSimilarityLinearNormalized( count, max_delta));
			  SCORES =  SCORES + weight * seq_sim;
		  }
		  else
		  {
			  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarityLinearNormalized( count, max_delta));

			  SCORES = weight * SCORES;
		  }
		  ++count;

      }
      else if( type.size() > 0)
      {
    	  std::cout << "WARNING: The type: " << type << "\n was not recognised!" << "\n";
      }
  }
  read.close();
  read.clear();
}
*/


void  ReadSimilarityScoreFileMsa
(
		ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORES,
		std::vector< Sequence> &MULTSEQALG1,
		std::vector< Sequence> &MULTSEQALG2,
		const std::string &SCORE_AND_PROFILE_FILE,
		std::set<char> &DEFINED_AMINO_ACIDS
)
{
  std::string

	line,
    file,
    type,
    first_file,
	second_file,
    tmp,
    file_ending,
    sliding_window_type;

	 std::vector< std::string>
		splitter;

  double
	  header_lines,
	  weight;

  size_t
	  column;

  int
	  splitter_size;

  size_t
	  count	=	0,
	  window_size;

  bool
	  file_not_empty = "FALSE";

  std::ifstream
	test_read,
    read;

  read.open( SCORE_AND_PROFILE_FILE.c_str());

	if( !read)
	{
		std::cerr << "ERROR: The score_and_profile_file [" << SCORE_AND_PROFILE_FILE << "] was not found";
		exit( -1);
	}


	while (read)
	{
		std::getline( read, line);
		if (line != "")
		{
			  file_not_empty = "TRUE";
			  splitter = mystr::SplitString(line, " ");
			  splitter_size = splitter.size();

#ifdef SECURE
			  if (splitter_size < 6)
			  {
				  std::cerr << "ERROR: A line in your similarity_score_file contains too less parameters! Please have a look at the manual on how to set up this file correctly!\n";
				  exit( -1);
			  }
#endif
			  weight = mystr::ConvertStringToNumericalValue< double>(splitter[1]);
#ifdef SECURE
			  if( !mystr::IsNumerical( splitter[1]) || (weight  < 0))
			  {
				  std::cerr << "ERROR: weight has to be a positive number, or zero, found: <" << splitter[1] << ">" << "\n";
				  exit( -1);
			  }
#endif
			  type = splitter[3];

			  if( type == "ScaleSimilarity")
			  {
#ifdef SECURE
				  if (splitter_size != 10)
				  {
					  std::cerr << "ERROR: In the similarity_score_file a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" file: <filename> min_cutoff <double> \n";
					  exit( -1);
				  }
#endif
				  file = splitter[5];
#ifdef SECURE
				  test_read.open(  file.c_str());

				  if( !test_read)
				  {
					  std::cerr << "ERROR: The file called   [" <<   file << "] for the type "<< type << " was not found";
					  exit( -1);
				  }
				  test_read.close();
				  test_read.clear();
#endif
				  sliding_window_type = splitter[7];


				  window_size = mystr::ConvertStringToNumericalValue<int>(splitter[9]);
#ifdef SECURE
				  if( !mystr::IsNumerical(splitter[9]) || ( window_size < 0) )
				  {
					  std::cerr<< "ERROR: windowsize value for UniveralProfileSimilarity has to be a positive number or zero, found: <"<< tmp << ">" << "\n";
					  exit( -1);
				  }
#endif
				  for( std::vector< Sequence>::iterator itr = MULTSEQALG1.begin(); itr != MULTSEQALG1.end(); ++itr)
				  {
#ifdef SECURE
			  			  if( itr->size() == 0)
			  			  {
			  				  std::cerr << "ERROR: It was not possible to read in a sequence (necessary for using scales) from one of the sequences you have submitted"<< "\n";
			  				  exit( -1);
			  			  }
			  			  if( itr->size() <= window_size)
			  			  {
			  				  std::cerr << "ERROR: Sequences have to be at least as long as the size of the window requested! \nLength of sequence <" << itr->size() << ">, length of window "<< window_size << "\n";
			  				  exit( -1);
			  			  }
#endif
			  			  ReadScaleIntoGappySequence( *itr, file, sliding_window_type, window_size, DEFINED_AMINO_ACIDS);
			  	 }
				  for( std::vector< Sequence>::iterator itr = MULTSEQALG2.begin(); itr != MULTSEQALG2.end(); ++itr)
				  {
#ifdef SECURE
			  			  if( itr->size() == 0)
			  			  {
			  				  std::cerr << "ERROR: It was not possible to read in a sequence (necessary for using scales) from one of the sequences you have submitted"<< "\n";
			  				  exit( -1);
			  			  }
			  			  if( itr->size() <= window_size)
			  			  {
			  				  std::cerr << "ERROR: Sequences have to be at least as long as the size of the window requested! \nLength of sequence <" << itr->size() << ">, length of window "<< window_size << "\n";
			  				  exit( -1);
			  			  }
#endif
			  			  ReadScaleIntoGappySequence( *itr, file, sliding_window_type, window_size, DEFINED_AMINO_ACIDS);
			  	 }
				  if( SCORES)
				  {
					  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
						  seq_sim( new ScoreProfileSimilarity( count));
					  SCORES =  SCORES + weight * seq_sim;
				  }
				  else
				  {
					  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarity( count));
					  SCORES = weight * SCORES;
				  }
				  ++count;
			  }
			  else if( type == "UniversalProfileSimilarity")
			  {
#ifdef SECURE
				  if (splitter_size != 10 )
				  {
					  std::cerr << "ERROR: In the similarity_score_file a line about " << type <<" does not fulfill the correct formatting rules. The format for  " << type <<" has to be: weight: <double> type: "<< type <<" file: <filename> min_cutoff <double> \n";
					  exit( -1);
				  }
#endif
				  column = mystr::ConvertStringToNumericalValue< size_t>(splitter[5]);
#ifdef SECURE
				  if( !mystr::IsNumerical( splitter[5]) || ( column  < 0) )
				  {
					  std::cerr<< "ERROR: column value for UniversalProfileSimilarity has to be a positive number or zero, found: <"<< tmp<< "> in " <<  SCORE_AND_PROFILE_FILE << "\n";
					  exit( -1);
				  }
#endif
				  header_lines = mystr::ConvertStringToNumericalValue< size_t>( splitter[7]);
#ifdef SECURE
				  if( !mystr::IsNumerical( splitter[7]) || ( header_lines < 0) )
				  {
					  std::cerr<< "ERROR: header_lines value for UniversalProfileSimilarity has to be a positive number or zero, found: <"<< tmp<< "> in " <<  SCORE_AND_PROFILE_FILE << "\n";
					  exit( -1);
				  }
#endif
				  file_ending = splitter[9];

				  if( file_ending.substr( 0, 1) != ".")
				  {
					  file_ending.insert( 0, ".");
				  }

				  for( std::vector< Sequence>::iterator itr = MULTSEQALG1.begin(); itr != MULTSEQALG1.end(); ++itr)
				  {
#ifdef SECURE
					  if( itr->size() == 0)
					  {
						  std::cerr << "ERROR: It was not possible to read in a sequence (necessary for using scales) from one of the sequences you have submitted"<< "\n";
						  exit( -1);
					  }
					  if( itr->size() <= window_size)
					  {
						  std::cerr << "ERROR: Sequences have to be at least as long as the size of the window requested! \nLength of sequence <" << itr->size() << ">, length of window "<< window_size << "\n";
						  exit( -1);
					  }
#endif
					  file = itr->GetFastaHeader() + file_ending;
					  ReadUniversalProfile( *itr, file,  column, header_lines );
				  }
				  for( std::vector< Sequence>::iterator itr = MULTSEQALG2.begin(); itr != MULTSEQALG2.end(); ++itr)
				  {
#ifdef SECURE
					  if( itr->size() == 0)
					  {
						  std::cerr << "ERROR: It was not possible to read in a sequence (necessary for using scales) from one of the sequences you have submitted"<< "\n";
						  exit( -1);
					  }
					  if( itr->size() <= window_size)
					  {
						  std::cerr << "ERROR: Sequences have to be at least as long as the size of the window requested! \nLength of sequence <" << itr->size() << ">, length of window "<< window_size << "\n";
						  exit( -1);
					  }
#endif
					  file = itr->GetFastaHeader() + file_ending;
					  ReadUniversalProfile( *itr, file,column, header_lines);
				  }
				  if( SCORES)
				  {
					  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
						  seq_sim( new ScoreProfileSimilarity( count));
					  SCORES =  SCORES + weight * seq_sim;
				  }
				  else
				  {
					  SCORES = ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >( new ScoreProfileSimilarity( count));
					  SCORES = weight * SCORES;
				  }
				  ++count;

			  }

			  else
			  {
				  std::cerr << "ERROR: Only ScaleSimilarity and UniversalProfileSimilarity are supported for the alignment of two averaged MSAs\n";
				  exit( -1);
			  }
		 }
	}
	if (!file_not_empty)
	{
		  std::cerr << "ERROR: Your similarity_score_file: <" << SCORE_AND_PROFILE_FILE << "> is empty or contains an empty line at the beginning \n";
		  exit( -1);
	}
}


// Procedure to read in fasta files!

std::string ReadSequence( Sequence &SEQUENCE, const std::string &FASTA_FILE, std::set<char> &DEFINED_AMINO_ACIDS)
{

  std::ifstream 
    read( FASTA_FILE.c_str());

  std::string
	  fasta_header,
	  line;

	if( !read)
	{
		std::cerr <<"ERROR: The fasta_file ["<< FASTA_FILE << "] was not found" << "\n";
		exit( -1);
	}

  std::getline( read, line);

	// reads the first line in the fasta-file and asserts if the line does not start with a "<"
	// the first line normally looks like: >PDBNAME:CAHINID|PDBID|CHAIN|SEQUENCE
  DebugWrite( "first letter fasta file: " << line.substr( 0, 1));

#ifdef SECURE
  if( line.substr( 0, 1) != ">")
  {
	  std::cerr << "ERROR: The fasta_file ["<< FASTA_FILE <<"] is not in the correct fasta format. ";
	  std::cerr << "A fasta file has to contain a header line starting with '>' followed by the sequence in the subsequent lines \n";
	  exit( -1);
  }
#endif

  fasta_header = line.substr( 1, line.size() - 1);
  SEQUENCE.SetFastaHeader( fasta_header);

  while( read)
  {
	  line.clear();
	  std::getline( read, line);
	  DebugWrite( line);
	  if( line.size() > 0)
	  {
#ifdef SECURE
		  if( line.substr( 0, 1) == ">")
		  {
			  std::cout << "ERROR: There should be only one sequence in the fasta_file  "<< FASTA_FILE << ". Another '>' was found! \n";
			  exit( -1);
		  }
#endif
			// goes through the string with all letters of a line
			// itr has the name of the amino-acid at a certain position
			// stores them in the aminco_acid_sequence (first or second)
			// GeneralizedAminoAcid( const char &TYPE) : m_Type( TYPE), scales()
			// scales will be empty in this case
		  for( std::string::const_iterator itr = line.begin(); itr != line.end(); ++itr)
		  {
			  if( *itr != ' ')
			  {
				  SEQUENCE.push_back( GeneralizedAminoAcid( *itr));

					// DEFINED_AMINO_ACIDS is a list in which all amino_acids that appear in the two sequences will be stored
					// therefore a check takes place if the actual amino acid is already in the list or not
					// if not than it will be added to the list
				  DEFINED_AMINO_ACIDS.insert(*itr);
			  }
		  }
	  }
  }

  read.close();
  read.clear();
  return fasta_header;
}

void ReadMultipleSequences( std::vector< Sequence> &SEQUENCES, const std::string &FASTA_FILE, std::set<char> &DEFINED_AMINO_ACIDS)
{
	std::vector< Sequence>::reverse_iterator r_itr;

  std::ifstream
    read( FASTA_FILE.c_str());

	if( !read)
	{
		std::cerr <<"ERROR: The msa_file ["<< FASTA_FILE << "] was not found" << "\n";
		exit( -1);
	}

  std::string
    line;

  while( read)
  {
	  line.clear();
	  std::getline( read, line);
	  if( line.size() == 0)
	  {
		  continue;
	  }

	  if( line.substr( 0, 1) == ">")
	  {
		  SEQUENCES.push_back( Sequence());
		  r_itr = SEQUENCES.rbegin();
		  r_itr->SetFastaHeader( line.substr( 1, line.size()-1));
	  }
	  else
	  {
#ifdef SECURE
		  if( SEQUENCES.size() == 0)
		  {
			  std::cerr << "ERROR: The msa_file ["<< FASTA_FILE <<"] is not in the correct fasta format. ";
			  std::cerr << "A fasta file has to contain a header line starting with '>' followed by the sequence in the subsequent lines \n";
			  exit( -1);
		  }
#endif
		  for( std::string::const_iterator itr = line.begin(); itr != line.end(); ++itr)
		  {

			  if( *itr != ' ')
			  {
				  r_itr->push_back( GeneralizedAminoAcid( *itr));
				  DEFINED_AMINO_ACIDS.insert(*itr);
			  }
		  }
	  }
    }
  read.close();
  read.clear();
}




