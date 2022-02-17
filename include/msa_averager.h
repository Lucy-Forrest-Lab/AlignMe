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
//!  AverageMsa
//! calculates an average profile from a multiple sequence alignment
//!
//!
//!
//!
//! @author: Kamil Khafizov, Rene Staritzbichler
//! @date: 1.4.2010
/////////////////////////////////////////////////////////////////////////


#ifndef AVERAGER_MSA_H_
#define AVERAGER_MSA_H_


inline
Sequence
AverageMsa
(
		const std::vector< Sequence> &SEQUENCES,
		const double &FRACTION_ALLOWED_GAPS,
		std::vector< int> &ERASED_POSITION,
		std::string &FILE
)
{
	size_t
		position = 0;

	size_t
		size_seqs = SEQUENCES[0].size(),
		profile_size = 0;

	// determine number of profiles by finding first non-gap GAA
	for( std::vector< GeneralizedAminoAcid>::const_iterator itr = SEQUENCES[0].begin(); itr != SEQUENCES[0].end(); ++itr)
	{
		if( itr->GetType() != '-')
		{
			profile_size = itr->GetProfiles().size();
			break;
		}
	}

#ifdef SECURE
	if( profile_size == 0)
	{
		std::cerr << "ERROR: All positions of the first sequence in your MSA from the file  < "<< FILE<< " > are skipped because they only represent columns of the MSA that are below the number of the fraction_of_allowed_gaps! \n";
		std::cerr << "However, the first sequence should be the representative sequence for your MSA. Please choose another sequence as being the representative one or change the value for the fraction_of_allowed_gaps \n";
		exit( -1);
	}
#endif
	Sequence
		sequence( size_seqs, GeneralizedAminoAcid( 'X', profile_size));

	std::vector< int>
		gaps( size_seqs, 0);

	Sequence::iterator
		result_itr;

	for( std::vector< Sequence>::const_iterator seqs_itr = SEQUENCES.begin(); seqs_itr != SEQUENCES.end(); ++seqs_itr)
	{
		position = 0;
		result_itr = sequence.begin();
		for( Sequence::const_iterator itr = seqs_itr->begin(); itr != seqs_itr->end(); ++itr, ++position, ++result_itr)
		{
			if( itr->GetType() == '-')
			{
				++gaps[ position];
			}
			else
			{
				result_itr->SumToProfiles( itr->GetProfiles());
			}
		}
	}


	Sequence
		result;
	// remove gappy generalized amino acids from sequence;
	int
		aa = 0,
		max_nr_gaps = int( FRACTION_ALLOWED_GAPS * SEQUENCES.size()),
		nr_non_gaps;

	std::vector< int>::const_iterator
		gap_itr = gaps.begin();

	for( Sequence::iterator itr = sequence.begin(); itr != sequence.end(); ++itr, ++gap_itr, ++aa)
	{
		if( *gap_itr <= max_nr_gaps)
		{
			nr_non_gaps = SEQUENCES.size() - *gap_itr;  // TODO: check
			itr->Normalize( double( nr_non_gaps));
			result.push_back( *itr);
		}
		else
		{
			ERASED_POSITION.push_back( result.size());
		}
	}

	//std::cout << "RESULT " << result << "\n";
	return result;
}



#endif /* AVERAGER_MSA_H_ */
