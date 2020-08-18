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

#include "sliding_window.h"
#include "string_functions.h"

// the itr goes through the amino-acid-sequence that has been read in (f.e. first_sequence or second_sequence)
// itr->GetType () returns the m_Type (which is in this case an AA) that is stored at the position to which itr points to
// This value will be searched in Scale_Map in which values corresponding to a certain AA are stored as pairs
// with -> second you have access to this corresponding value
// This value will than be added to the sequence hat the position that itr points to
// Because of AddNewProfile(value) the value will be pushed back in m_Profiles of the corresponding GAA

void ReadScaleIntoSequence(AASequence &SEQUENCE,const std::map<char, double> &SCALE_MAP)
{
	DebugWrite(__FUNCTION__);
	for(  std::vector< GeneralizedAminoAcid>::iterator itr = SEQUENCE.begin(); itr != SEQUENCE.end(); ++itr)
    {
		// std::map<>::find() returns an iterator to key/value std::pair // has to be used rather then operator [] because map is passed as const
		itr->AddNewProfile( SCALE_MAP.find( itr->GetType())->second);
		// (*itr).AddNewProfile( SCALE_MAP.find( (*itr).GetType()->second);
    }
}


//  for(  std::vector< GeneralizedAminoAcid>::iterator itr = SEQUENCE.begin(); itr != SEQUENCE.end(); ++itr)
//  {
//	  GeneralizedAminoAcid gaa( *itr);
//	  char type( gaa.GetType());
//	  double scale_value( SCALE_MAP.find( type));
//	  gaa.AddNewProfile( value);
//    *itr = gaa;
//  }

// called from ReadProfilesAndScales
// with ReadScale( FIRST, file, sliding_window_type, window_size);
// with ReadScale( SECOND, file, sliding_window_type, window_size);
// file is the scale_file, f.e. KD.txt
void ReadScale
(
		AASequence &SEQUENCE,
		const std::string &SCALE_FILE,
		const std::string &SLIDING_WINDOW_TYPE,
		const int &WINDOW_SIZE,
		std::set<char> &DEFINED_AMINO_ACIDS
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
				std::cout << "WARNING: Amino acid <" << aa_type << "> appears multiple times in the scale (" << SCALE_FILE << "). Its previous value: <" << scale[ aa_type] << ">, has been overwritten and is now: <" << value_str << "\n";
			}
			scale.insert( std::make_pair( aa_type, scale_value));
		}
		else
		{
			 std::cout << "WARNING: the value "<< aa_type <<" from the scale " << SCALE_FILE << " does not appear in any of your sequences " << "\n";
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
		AASequence &SEQUENCE,
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
			std::cerr << "ERROR: Only single-letter amino acid keys are allowed" << "\n";
			std::cerr << "ERROR: found " << aa_str << " in " << SCALE_FILE << "\n";
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
				std::cout << "WARNING: amino acid <" << aa_type << "> appears multiple times in (" << SCALE_FILE << "), previous value: <" << scale[ aa_type] << ">, is now: <" << value_str << "\n";
			}
			scale.insert( std::make_pair( aa_type, scale_value));
		}
		else
		{
			std::cout << "WARNING: the value "<< aa_type <<" from the scale " << SCALE_FILE << " does not appear in any of your sequences " << "\n";
		}
	}
	// check if the amino_acids that are in the sequences that have been read in, are in the scale
	// if an amino acid is in the sequence and not in the scale, than the program will quit
	for (std::set<char>::iterator acid_itr = DEFINED_AMINO_ACIDS.begin(); acid_itr!= DEFINED_AMINO_ACIDS.end(); ++acid_itr)
	{
		std::map<char, double>::const_iterator itr = scale.find(*acid_itr);
		if( *acid_itr != '-' && itr == scale.end())
		{
			std::cerr << "ERROR: the amino acid " << *acid_itr<< " does not exist in the scale " << SCALE_FILE<< " that you want to use" << "\n";
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
	  std::cerr << "ERROR: for multiple sequence alignments the only window type that is currently allowed is 'triangular_msa'" << "\n";
	  exit( -1);
  }
  read.close();
  read.clear();
}


// Is called by ReadProfile( FIRST, first_file); or ReadProfile( SECOND, first_file)
// where FIRST and SECOND are constructors of GeneralizedAminoAcids
void  ReadProfile( AASequence &SEQUENCE, const std::string &PROFILE_FILE)
{
//	std::cout <<  __FUNCTION__ << "\n";
	std::ifstream
		read( PROFILE_FILE.c_str());
	double
		value;
	std::string
		str;

	if( !read)
	{
		if ( PROFILE_FILE != "")
		{
			std::cerr << "ERROR: The profile file [ " << PROFILE_FILE << "] was not found" << "\n";
		}
		else
		{
			std::cerr << "ERROR: You have to provide 2 profile files to make an alignment based on ProfileSimilarity" << "\n";
		}
		exit( -1);
	}

	// If the GeneralizedAminoAcidSequence is empty an nothing has been added to it yet
	// then the profile_file will be opened and as long as there are values been read in:
	// m_Type (X) will be added because of GeneralizedAminoAcid()
	// m_Profiles (value) will also be added at this position
	if( SEQUENCE.size() == 0)
	{
//		std::cout <<  __FUNCTION__ << " new sequence" << "\n";
		while( read)
		{
#ifdef SECURE
			str.clear();
			read >> str;

			if( str.size() == 0)
			{
				continue;
			}

			if( !mystr::IsNumerical( str))
			{
				std::cerr << "ERROR: The entry " << str << " is not numerical" << "\n";
				exit( -1);
			}

			value = mystr::ConvertStringToNumericalValue< double>( str);

			SEQUENCE.push_back( GeneralizedAminoAcid());
			SEQUENCE.back().AddNewProfile( value);
#else
			value = std::numeric_limits< double>::max();
			read >> value;
			if( value != std::numeric_limits< double>::max())
			{
				SEQUENCE.push_back( GeneralizedAminoAcid());
				SEQUENCE.back().AddNewProfile( value);
			}
#endif
		}
	}

	// If there are already values stored in the GAA-Sequence:
	// the itr goes through the whole sequence and points to a certain GAA
	// a value corresponding to the profile_file will be added to this certain GAA in m_Profiles
	// m_Profiles can contain a list of different values corresponding to a certain GAA
	else
	{
//		std::cout <<  __FUNCTION__ << " sequence exists" << "\n";
		std::vector< GeneralizedAminoAcid>::iterator itr = SEQUENCE.begin();
		while( !read.eof() && itr != SEQUENCE.end())
		{
#ifdef SECURE
			str.clear();
			read >> str;
			if( str.size() > 0)
			{
				if( !mystr::IsNumerical( str))
				{
					std::cerr << "ERROR: The entry " << str << " is not numerical" << "\n";
					exit( -1);
				}
				itr->AddNewProfile(mystr::ConvertStringToNumericalValue<double>(str));
			}
#else
			str.clear();
			read >> str;
			itr->AddNewProfile(mystr::ConvertStringToNumericalValue<double>(str));
#endif

		  ++itr;
		}

		read >> str;

		if( !read.eof())
		{
			std::cout << "ERROR: The profile file [" << PROFILE_FILE << "] contains less values than there are amino acids in the corresponding sequence that you submitted!" << "\n";
			std::cout << "ERROR: Please be sure that you have submitted the correct files and colum-number" << "\n";
			exit(-1);
		}

		if( itr != SEQUENCE.end())
		{
			std::cout << "ERROR: The profile file [" << PROFILE_FILE << "] contains more values than there are amino acids in the corresponding sequence that you submitted!" << "\n";
			std::cout << "ERROR: Please be sure that you have submitted the correct files and column number" << "\n";
			exit(-1);
		}

	}
	read.close();
	read.clear();
}


// Is called by ReadProfile( FIRST, first_file); or ReadProfile( SECOND, first_file)
// where FIRST and SECOND are constructors of GeneralizedAminoAcids
void  ReadUniversalProfile( AASequence &SEQUENCE, const std::string &PROFILE_FILE, const int &COLUMN, const int  &HEADER_LINES)
{
//	std::cout <<  __FUNCTION__ << "\n";
	std::ifstream
		read( PROFILE_FILE.c_str());

	std::string
		tmp,
		line;

	double
		value;

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
		std::cout <<  "headerlines are " << HEADER_LINES << std::endl;
		std::cout <<  "skipped lines are " << line << std::endl;
	}
	// If the GeneralizedAminoAcidSequence is empty an nothing has been added to it yet
	// then the profile_file will be opened and as long as there are values been read in:
	// m_Type (X) will be added because of GeneralizedAminoAcid()
	// m_Profiles (value) will also be added at this position
	if (SEQUENCE.size() == 0)
	{
//		std::cout <<  __FUNCTION__ << " new sequence" << "\n";
		while (read)
		{
			value = std::numeric_limits<double>::max();
			tmp = "";

			for (int i = 0; i != COLUMN; i++)
			{
				read >> tmp;
			}

			if( tmp != "")
			{
#ifdef SECURE
				if( !mystr::IsNumerical( tmp))
				{
					std::cerr << "ERROR: The value: <" << tmp << "> in column " << COLUMN << " of profile file (" << PROFILE_FILE << ") is not numerical" << "\n";
					exit( -1);
				}
#endif
				value = mystr::ConvertStringToNumericalValue< double>( tmp);
				SEQUENCE.push_back(GeneralizedAminoAcid());
				SEQUENCE.back().AddNewProfile(value);
			}
		   std::getline( read, line);
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
			tmp = "";

			for (int i = 0; i != COLUMN; i++)
			{
				read >> tmp;
			}

			if( tmp != "")
			{
#ifdef SECURE
				if( !mystr::IsNumerical( tmp))
				{
					std::cerr << "ERROR: The value: <" << tmp << "> in column " << COLUMN << " of profile file (" << PROFILE_FILE << " is not numerical" << "\n";
					exit( -1);
				}
#endif
				value = mystr::ConvertStringToNumericalValue<double>( tmp);
				itr->AddNewProfile( value);
				++itr;
				std::getline( read, line);
			}
		}

		read >> tmp;

#ifdef SECURE
		if( !read.eof())
		{
			std::cerr << "ERROR: The profile in " << PROFILE_FILE << " is longer than the corresponding sequence that you submitted!" << "\n";
			std::cerr << "ERROR: If you submit sequences and profiles at the same time, please be sure that they correspond to one another" << "\n";
			exit( -1);
		}
		if( itr != SEQUENCE.end())
		{
			std::cerr << "ERROR: The profile in " << PROFILE_FILE << " is shorter than the corresponding sequence that you submitted!" << "\n";
			std::cerr << "ERROR: If you submit sequences and profiles at the same time, please be sure that they correspond to one another" << "\n";
			exit( -1);
		}
#endif
	}
	read.close();
	read.clear();
}




// Called from sweet_aling.cpp with   ReadProfilesAndScales( first_sequence, second_sequence, score_and_profile_file);
// Score_and_Profile-file looks like: 10  ScaleSimilarity  KD.txt triangular 13 and can have several command-lines
// Calculates a score for each position

void  ReadProfilesAndScales
(
		AASequence &FIRST,
		AASequence &SECOND,
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
    file,
    type,
    first_file,
    second_file,
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
	  if (tmp == "none" || tmp == "")
	  {
		  continue;

	  }
	  if( !mystr::IsNumerical( tmp) || (weight = mystr::ConvertStringToNumericalValue<int>( tmp) < 0))
	  {
		  std::cerr << "ERROR: weight has to be a positive number, or zero, found: <" << tmp << ">" << "\n";
		  exit( -1);
	  }
#else
	  read >> tmp >> weight;
#endif

	  read >> tmp >> type;

	  //  weight similaritytype (f.e. BLOSUM or PAM)
	  if( type == "SequenceSimilarity")
	  {
		  read >> tmp >> file;
	  }
		// weight similaritytype firstprofile_file secondprofile_file
		// FIRST = first sequence  SECOND = second sequence
	  else if( type == "ProfileSimilarity") //TODO: get rid of it later as it duplicates UniverslaProfileSimilarity
	  {
		  read >> tmp >> first_file;
		  ReadProfile( FIRST, first_file);
		  read >> tmp >> second_file;
		  ReadProfile( SECOND, second_file);
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
		  read >> tmp >> first_file;
		  ReadUniversalProfile( FIRST, first_file, column, header_lines);
		  read >> tmp >> second_file;
		  ReadUniversalProfile( SECOND, second_file, column, header_lines);

		  WRITE_PROFILES_HEADER1.push_back(" based on profile " + first_file + ", extracted column: "+mystr::NumericalValueToString(column)+"\n");
		  WRITE_PROFILES_HEADER2.push_back(" based on profile " + second_file + ", extracted column: "+mystr::NumericalValueToString(column)+"\n");

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
#ifdef SECURE
		  if (sliding_window_type == "")
		  {
			  std::cerr << "ERROR: You have to provide a windowtype and a windowsize if you want to do an alignment based on scale similarity"<<"\n";
			  exit(-1);
		  }
#endif
#ifdef SECURE
		  if( FIRST.size() == 0 || SECOND.size() == 0)
		  {
			  std::cerr << "ERROR: You have to provide two sequences if you want to make an alignment based on scales"<<"\n";
			  std::cerr << "ERROR: Length of provided sequences, first sequence: <" << FIRST.size() << ">, second sequence: <" << SECOND.size() << ">" << "\n";
			  exit( -1);
		  }
		  tmp.clear();
		  read >> tmp >> tmp;
		  if (tmp == "")
		  {
			  std::cerr << "ERROR: You have to provide a windowsize if you want to do an alignment based on scale similarity"<<"\n";
			  exit(-1);
		  }
		  if( !mystr::IsNumerical( tmp) || ( window_size = mystr::ConvertStringToNumericalValue<int>(tmp)) < 0)
		  {
			  std::cerr<< "ERROR: windowsize value for UniveralProfileSimilarity has to be a positive number or zero, found: <"<< tmp << ">" << "\n";
			  exit( -1);
		  }
		  if( int(FIRST.size()) <= window_size || int(SECOND.size()) <= window_size)
		  {
			  std::cerr << "ERROR: Sequences have to be at least as long as the size of the window requested! "<< "\n";			  std::cerr << "ERROR: Length of provided sequences, first: <" << FIRST.size() << ">, second: <" << SECOND.size() << ">" << "\n";
			  std::cerr << "ERROR: Length of provided sequences, first: <" << FIRST.size() << ">, second: <" << SECOND.size() << ">" << "\n";
			  exit( -1);
		  }
#else
		  read >> tmp >> window_size;
#endif
		  ReadScale( FIRST, file, sliding_window_type, window_size, DEFINED_AMINO_ACIDS);
		  ReadScale( SECOND, file, sliding_window_type, window_size, DEFINED_AMINO_ACIDS);

		  WRITE_PROFILES_HEADER1.push_back(" based on scale " +  file + " for sequence "+  FIRST_FASTA_ID + " using windowtype " + sliding_window_type + " with size " + mystr::NumericalValueToString(window_size)+"\n");
		  WRITE_PROFILES_HEADER2.push_back(" based on scale " +  file + " for sequence "+  SECOND_FASTA_ID + " using windowtype "+ sliding_window_type + " with size " + mystr::NumericalValueToString(window_size)+"\n");
	  }
	  else if( type.size() > 0)
	  {
		  std::cerr << "ERROR: similarity type requested is not recognised:<" << type << ">" << "\n";
		  exit( -1);
	  }
  }
  read.close();
  read.clear();
}




void  ReadProfilesAndScalesMsa
(
		AASequence &SEQ,
		const std::string &SCORE_AND_PROFILE_FILE,
		std::set<char> &DEFINED_AMINO_ACIDS
)
{
  std::string
    file,
    type,
    first_file,
    tmp,
    sliding_window_type;

  double
    weight;

  size_t
	window_size;

  std::ifstream
    read( SCORE_AND_PROFILE_FILE.c_str());

	if( !read)
	{
		std::cerr << "ERROR: The score_and_profile_file [" <<SCORE_AND_PROFILE_FILE << "] was not found";
		exit( -1);
	}

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
	  if( !mystr::IsNumerical( tmp) || ( weight = mystr::ConvertStringToNumericalValue< double>( tmp)) < 0)
	  {
		  std::cerr << "ERROR: weight has to be positive number or zero, found: <" << tmp << ">" << "\n";
		  exit( -1);
	  }
#else
	  read >> tmp >> weight;
#endif
	  read >> tmp >> type;

	  if( type == "ScaleSimilarity")
	  {
		  read >> tmp >> file;
		  read >> tmp >> sliding_window_type;

#ifdef SECURE
		  if( SEQ.size() == 0)
		  {
			  std::cerr << "ERROR: It was not possible to read in a sequence (necessary for using scales) from one of the sequences you have submitted"<< "\n";
			  exit( -1);
		  }
		  read >> tmp >> tmp;
		  if( !mystr::IsNumerical( tmp) || ( window_size = mystr::ConvertStringToNumericalValue< size_t>( tmp)) < 0)
		  {
			  std::cerr<< "ERROR: header_lines value for UniveralProfileSimilarity has to be a positive number or zero, found: <"<< tmp << ">" << "\n";
			  exit( -1);
		  }
		  if( SEQ.size() <= window_size)
		  {
			  std::cerr << "ERROR: Sequences have to be at least as long as the size of the window requested! \nLength of sequence <" << SEQ.size() << ">, length of window "<< window_size << "\n";
			  exit( -1);
		  }
#else
		  read >> tmp >> window_size;
#endif

		  ReadScaleIntoGappySequence( SEQ, file, sliding_window_type, window_size, DEFINED_AMINO_ACIDS);

	  }
	  else if( type == "SequenceSimilarity" || type == "ProfileSimilarity" || type == "UniversalProfileSimilarity")
	  {
		  std::cerr << "ERROR: Only 'ScaleSimilarity' is defined at the moment to align multiple sequence alignments" << "\n";
		  exit(-1);
	  }
  }
  read.close();
  read.clear();
}
//////////////////////////////////////////////////////////////////////

// is called by ReadScores( scores, score_and_profile_file) from sweet_align.cpp
// Score_and_Profile-file looks like: 10  ScaleSimilarity  KD.txt triangular 13 and can have several command-lines

void ReadScores
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
	  if(tmp == "none" || tmp == "")
	  {
		  continue;
	  }
	  if( tmp != "" && ( !mystr::IsNumerical( tmp) || ( weight = mystr::ConvertStringToNumericalValue< double>( tmp)) < 0))
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
    	  std::cout << "WARNING: The type: " << type << "\n was not recognised!" << "\n";
      }
  }
  read.close();
  read.clear();
}



std::string ReadSequence( AASequence &SEQUENCE, const std::string &FASTA_FILE, std::set<char> &DEFINED_AMINO_ACIDS)
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
  if( line.substr( 0, 1) != ">")
  {
	  std::cerr << "ERROR: The file ["<< FASTA_FILE <<"] is not a correct fasta file" << "\n";
	  std::cerr << "ERROR: A fasta file has to contain a sequence, and a header line starting with '>'" << "\n";
	  exit( -1);
  }

  fasta_header = line.substr( 1, line.size() - 1);

  while( read)
  {
	  line.clear();
	  std::getline( read, line);
	  DebugWrite( line);
	  if( line.size() > 0)
	  {
		  if( line.substr( 0, 1) == ">")
		  {
			  std::cout << "WARNING: There should be only one sequence in "<< FASTA_FILE << ". (another '>' was found)" << "\n";
			  std::cout << "WARNING: Only the first sequence of "<< FASTA_FILE << " will be used" << "\n";
			  return "";
		  }
			// goes through the string with all letters of a line
			// itr has the name of the amino-acid at a certain position
			// stores them in the aminco_acid_sequence (first or second)
			// GeneralizedAminoAcid( const char &TYPE) : m_Type( TYPE), scales()
			// scales will be empty in this case
		  for( std::string::const_iterator itr = line.begin(); itr != line.end(); ++itr)
		  {
			  if( *itr != ' ')
			  {
//				  DebugWrite( *itr);
				  SEQUENCE.push_back( GeneralizedAminoAcid( *itr));

					// DEFINED_AMINO_ACIDS is a list in which all amino_acids that appear in the two sequences will be stored
					// therefore a check takes place if the actual amino acid is already in the list or not
					// if not than it will be added to the list
				  DEFINED_AMINO_ACIDS.insert(*itr);
			  }
		  }
		  DebugWrite( "  <========");
	  }
  }

  read.close();
  read.clear();
  return fasta_header;
}


void ReadMultipleSequences( std::vector< AASequence> &SEQUENCES, const std::string &FASTA_FILE, std::set<char> &DEFINED_AMINO_ACIDS, bool EQUILENGTH = true)
{
	std::vector< AASequence>::reverse_iterator r_itr;

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
		  SEQUENCES.push_back( AASequence());
		  r_itr = SEQUENCES.rbegin();
	  }
	  else
	  {
		  if( SEQUENCES.size() == 0)
		  {
			std::cerr << "ERROR: The msa_file ["<< FASTA_FILE << "] has to be in fasta-format starting with a > followed by a header and a amino acid sequece in the next lines" << "\n";
			exit( -1);
		  }
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

#ifdef SECURE
  if( EQUILENGTH)
  {
	int
		i = 0,
		length = SEQUENCES[0].size();
	for( std::vector< std::vector< GeneralizedAminoAcid> >::const_iterator itr = SEQUENCES.begin() + 1; itr != SEQUENCES.end(); ++itr, ++i)
	{
		if( (int) itr->size() != length)
		{
			std::cerr << "ERROR: Your multiple sequence alignment contains sequences of different length \n ";
			std::cerr << "ERROR: Length of sequence "<< i+2 << " is " << itr->size() << " is not matching the reference length " << length << " given by the first sequence" << "\n";
			std::cerr << "ERROR: All sequences within the MSA need to have the same length!\n";
			exit( -1);
		}
	}
  }
#endif

  read.close();
  read.clear();
}

#endif /* READER_H */

