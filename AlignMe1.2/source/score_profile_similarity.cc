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
//!  A collection of functions to evaluate the similarity of two profiles.
//! The ScoreProfileSimilarity class returns a sum over the differences.
//!
//! TODO: seperate functions from class
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////

#include "../include/score_profile_similarity.h"

#include <cmath>
#include "../include/macro_functions_read_write.h"



ScoreProfileSimilarity::ScoreProfileSimilarity()
: m_IdInProfile( std::numeric_limits< size_t>::max())
{}

ScoreProfileSimilarity::ScoreProfileSimilarity( const size_t &ID)
: m_IdInProfile( ID)
{
  DebugWrite( __FUNCTION__);
}

ScoreProfileSimilarity::~ScoreProfileSimilarity(){}


double ScoreProfileSimilarity::operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
{
//	std::cout << "id: " << m_IdInProfile << std::endl;
//	std::cout << AA.first.GetProfiles()[ m_IdInProfile] << std::endl;
//	std::cout << AA.second.GetProfiles().size() << std::endl;
	double d = -fabs( AA.first.GetProfiles()[ m_IdInProfile] - AA.second.GetProfiles()[ m_IdInProfile]);
	//std::cout << d << "_";
	//   DebugWrite( __FUNCTION__ << "  " << d);
	return d;
}

std::ostream &ScoreProfileSimilarity::Write( std::ostream &STREAM) const
{
  STREAM << "ScoreProfileSimilarity::Write()" << "\n";
  STREAM << "profile-id: " << m_IdInProfile << "\n";
  return STREAM;
}


int NrAligned( const std::vector< std::pair< int, int> > &ALIGNMENT)
{
	int count = 0;
	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max())
		{
			++count;
		}
	}
	return count;
}

int LengthAligned( const std::vector< std::pair< int, int> > &ALIGNMENT)
{
	int length = 0;
	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		++length;
	}
	return length;
}

// VERSION 1.1================================ till end of file !
//This score is good when aligning sequences of similar length and when profiles match well. Do not care explicitly about gaps.
double PdsScoreAlignmentLength( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND, const double &RATIO)
{
	assert( RATIO >= 0.0 && RATIO <= 1.0);
	double score( 0.0);
	int aligned_count( 0);
	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max())
		{
			std::vector< double>
				first_profile( FIRST[ itr->first].GetProfiles()),
				second_profile( SECOND[ itr->second].GetProfiles());
			for( std::vector< double>::const_iterator first_itr = first_profile.begin(), second_itr = second_profile.begin(); first_itr != first_profile.end() && second_itr != second_profile.end(); ++first_itr, ++second_itr)
			{
				score += /* weight * */ pow( *first_itr - *second_itr, 2.0); // TODO: introduce weights!
			}
			++aligned_count;
		}
	}
	if( FIRST.size() > 0 && SECOND.size() > 0)
	{
		if( aligned_count >= RATIO * ALIGNMENT.size())
		{
			return  sqrt( score) * ALIGNMENT.size() / aligned_count ;
		}
		else
		{
			return 55555.55555;
		}
	}
	return std::numeric_limits< double>::max();
}

//This score is good when profiles match well and do not care of differences in length. Do not care explicitly about gaps.
double PdsScoreAligned( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND, const double &RATIO)
{
	assert( RATIO >= 0.0 && RATIO <= 1.0);
	double score( 0.0);
	int aligned_count( 0);
	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max())
		{
			std::vector< double>
				first_profile( FIRST[ itr->first].GetProfiles()),
				second_profile( SECOND[ itr->second].GetProfiles());
			for( std::vector< double>::const_iterator first_itr = first_profile.begin(), second_itr = second_profile.begin(); first_itr != first_profile.end() && second_itr != second_profile.end(); ++first_itr, ++second_itr)
			{
				score += /* weight * */ pow( *first_itr - *second_itr, 2.0); // TODO: introduce weights!
			}
			++aligned_count;
		}
	}
	if( FIRST.size() > 0 && SECOND.size() > 0)
	{
		if( aligned_count >= RATIO * ALIGNMENT.size())
		{
			return  sqrt( score) / aligned_count;
		}
		else
		{
			return 55555.55555;
		}
	}
	return std::numeric_limits< double>::max();
}


//This score is good when profiles match well and few gaps were found.
double PdsScoreAlignedAndGaps( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND, const double &RATIO)
{
	assert( RATIO >= 0.0 && RATIO <= 1.0);
	double score( 0.0);
	int aligned_count( 0);
	int gap_count( 0);

	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max())
		{
			std::vector< double>
				first_profile( FIRST[ itr->first].GetProfiles()),
				second_profile( SECOND[ itr->second].GetProfiles());
			for( std::vector< double>::const_iterator first_itr = first_profile.begin(), second_itr = second_profile.begin(); first_itr != first_profile.end() && second_itr != second_profile.end(); ++first_itr, ++second_itr)
			{
				score += /* weight * */ pow( *first_itr - *second_itr, 2.0); // TODO: introduce weights!
			}
			++aligned_count;
		}
		else
		{
			++gap_count;
		}
	}
	if( FIRST.size() > 0 && SECOND.size() > 0)
	{
		if( aligned_count >= RATIO * ALIGNMENT.size())
		{
			return  sqrt( score) * (gap_count + 1) / (aligned_count + 1); // + 1 is to avoid problems when gap_count
		}
		else
		{
			return 55555.55555;
		}
	}
	return std::numeric_limits< double>::max();
}


double PdsScoreAlignedMin( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND)
{
	double score( 0.0);
	int aligned_count( 0);
	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max())
		{
			std::vector< double>
				first_profile( FIRST[ itr->first].GetProfiles()),
				second_profile( SECOND[ itr->second].GetProfiles());
			for( std::vector< double>::const_iterator first_itr = first_profile.begin(), second_itr = second_profile.begin(); first_itr != first_profile.end() && second_itr != second_profile.end(); ++first_itr, ++second_itr)
			{
				score += /* weight * */ pow( *first_itr - *second_itr, 2.0); // TODO: introduce weights!
			}
			++aligned_count;
		}
	}
	if( FIRST.size() > 0 && SECOND.size() > 0)
	{
		return  sqrt( score) * std::min( FIRST.size(), SECOND.size()) / aligned_count ;
	}
	return std::numeric_limits< double>::max();
}

double PdsScoreAlignedMax( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND)
{
	double score( 0.0);
	int aligned_count( 0);
	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max())
		{
			std::vector< double>
				first_profile( FIRST[ itr->first].GetProfiles()),
				second_profile( SECOND[ itr->second].GetProfiles());
			for( std::vector< double>::const_iterator first_itr = first_profile.begin(), second_itr = second_profile.begin(); first_itr != first_profile.end() && second_itr != second_profile.end(); ++first_itr, ++second_itr)
			{
				score += /* weight * */ pow( *first_itr - *second_itr, 2.0); // TODO: introduce weights!
			}
			++aligned_count;
		}
	}
	if( FIRST.size() > 0 && SECOND.size() > 0)
	{
		return  sqrt( score) * std::max( FIRST.size(), SECOND.size()) / aligned_count ;
	}
	return std::numeric_limits< double>::max();
}

//This score is good when aligning sequences of similar length (because of gaps!) and when profiles match well. Mean value for gaps.
double PdsScore( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND, const double &GAP_EXTENSION_PENALTY, const double &RATIO)
{
	double score( 0.0);
	int aligned_count( 0);
	int gap_count (0);
	double value( 0.0);
	double mean( 0.0);
	std::vector< double> differences;

	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first == std::numeric_limits< int>::max() || itr->second == std::numeric_limits< int>::max())
		{
			gap_count++;
			//score += pow( GAP_EXTENSION_PENALTY, 2.0); //TODO: use this??? square this?
		}
		else
		{
			std::vector< double>
				first_profile( FIRST[ itr->first].GetProfiles()),
				second_profile( SECOND[ itr->second].GetProfiles());
				for( std::vector< double>::const_iterator first_itr = first_profile.begin(), second_itr = second_profile.begin(); first_itr != first_profile.end() && second_itr != second_profile.end(); ++first_itr, ++second_itr)
			{
					value = *first_itr - *second_itr;
					score += /* weight * */ pow( value, 2.0); // TODO: introduce weights!
					differences.push_back( fabs( value));

			}
			++aligned_count;
		}
	}

	for( std::vector< double>::const_iterator itr = differences.begin(); itr != differences.end(); ++itr)
	{
		mean += *itr;
	}
	mean /= differences.size();

	score += gap_count * pow( mean, 2.0);

	if( aligned_count >= RATIO * ALIGNMENT.size())
	{
		return sqrt( score / aligned_count);
	}
	return 77777.77777 ;//std::numeric_limits< double>::max();
}


double PdsFixedGaps( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND, const double &GAP_EXTENSION_PENALTY)
{
	double score( 0.0);
	int aligned_count( 0);
	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first == std::numeric_limits< int>::max() || itr->second == std::numeric_limits< int>::max())
		{
			score += pow( GAP_EXTENSION_PENALTY, 2.0); //TODO: use this??? square this?
		}
		else
		{
			std::vector< double>
				first_profile( FIRST[ itr->first].GetProfiles()),
				second_profile( SECOND[ itr->second].GetProfiles());
				for( std::vector< double>::const_iterator first_itr = first_profile.begin(), second_itr = second_profile.begin(); first_itr != first_profile.end() && second_itr != second_profile.end(); ++first_itr, ++second_itr)
			{
				score += /* weight * */ pow( *first_itr - *second_itr, 2.0); // TODO: introduce weights!
			}
			++aligned_count;
		}
	}

	if( aligned_count != 0)
	{
		return sqrt( score / aligned_count);
	}
	return std::numeric_limits< double>::max();
}





std::vector< double> CalculatePearsonCorrelations( const std::vector< std::pair< int, int> > &ALIGNMENT, const Sequence &FIRST, const Sequence &SECOND, const double &RATIO)
{
	std::vector< double>
		result;
	double
		sum_x( 0.0),
		sum_y( 0.0),
		sum_x_square( 0.0),
		sum_y_square( 0.0),
		match_count( 0.0),
		sum_x_times_y( 0.0),
		shortest_seq_length( 0.0),
		length_first( FIRST.size()),
		length_second( SECOND.size()),
		length( ALIGNMENT.size()),
		gap_count( 0.0);

	shortest_seq_length = length_first < length_second ?  length_first : length_second;

	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( ( itr->first == std::numeric_limits< int>::max() || itr->second == std::numeric_limits< int>::max()))
		{
			++gap_count;
		}
		else
		{

			std::vector< double>
			  first_profile( FIRST[ itr->first].GetProfiles()),
			  second_profile( SECOND[ itr->second].GetProfiles());

			std::vector< double>::const_iterator first_itr = first_profile.begin();

			for( std::vector< double>::const_iterator first_itr = first_profile.begin(), second_itr = second_profile.begin(); first_itr != first_profile.end() && second_itr != second_profile.end(); ++first_itr, ++second_itr)
			{
				sum_x += *first_itr;
				sum_y += *second_itr;
				sum_x_square += pow( *first_itr, 2.0);
				sum_y_square += pow( *second_itr, 2.0);
				sum_x_times_y += *first_itr * *second_itr;
			}
			++match_count;
		}
	}
	sum_x = sum_x/3;
	sum_y  = sum_y/3;
	sum_x_square = sum_x_square/3;
	sum_y_square  = sum_y_square/3;
	sum_x_times_y = sum_x_times_y/3;

	double correlation( ( length * sum_x_times_y - ( sum_x * sum_y))
		/ ( sqrt( length * sum_x_square - pow( sum_x, 2.0)) * sqrt( length * sum_y_square - pow( sum_y, 2))));

	if( match_count < RATIO * shortest_seq_length)
	{
		correlation = 1111.1111;
	}
	if( gap_count == 0.0)
	{
		gap_count = 1e-4;
	}

	result.push_back( correlation);
	result.push_back( match_count * correlation);
	result.push_back( match_count / shortest_seq_length * correlation);
	result.push_back( match_count / double( ALIGNMENT.size()) * correlation);
	result.push_back( match_count / gap_count * correlation);

	return result;
}

