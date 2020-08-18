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
std::pair< double, std::vector< std::pair< int, int> > >
ExtractPairAlignmentFromAlignedMSA
(
	const std::pair< double, std::vector< std::pair< int, int> > >  &ALIGNED_AMSA,
	const Sequence &FIRST_MSA_SEQ,
	const std::vector< int> &FIRST_IGNORED_POSITIONS,
	const Sequence &SECOND_MSA_SEQ,
	const std::vector< int> &SECOND_IGNORED_POSITIONS
)
{
	std::vector< std::pair< int, int> >
		aligned_list;
	Sequence
		first_sequence = FIRST_MSA_SEQ,  // selected sequence with window averaged scale values
		second_sequence = SECOND_MSA_SEQ;

#ifdef DEBUG
	std::ofstream write( "extracted_sequence_from_first_MSA.amf");
	write << first_sequence;
	write.close();

	write.open( "extracted_sequence_from_second_MSA.amf");
	write << second_sequence;
	write.close();
#endif

	int
		tmp,
		msa_a_id,
		msa_b_id,
		local_first_shift,
		local_second_shift,
		first_shift = 0,
		second_shift = 0,
		previous_first_shift = 0,
		previous_second_shift = 0,
		max = std::numeric_limits< int>::max();

	std::vector< int>
		first_ignored_positions = FIRST_IGNORED_POSITIONS,
		second_ignored_positions = SECOND_IGNORED_POSITIONS;

	// loop through the alignment of the 2 averaged MSAs, adding missing information about the query sequences
	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNED_AMSA.second.begin(); itr != ALIGNED_AMSA.second.end(); ++itr)
	{
		// indices of the averaged MSAs
		msa_a_id = itr->first;
		msa_b_id = itr->second;

//		std::cout << "now: pos " << /*cc++ <<*/ ": " << msa_a_id << " - " << msa_b_id << "\n";

		// adjust shift here

		local_first_shift = 0;
		local_second_shift = 0;

		if( msa_a_id < max) // max is gap value (gap created in the Alignment of AveragedMSA)
		{
			while
			(
					first_ignored_positions.size() > 0 &&      // check that the vector is not yet shrunk to zero
					first_ignored_positions[0] <= msa_a_id     // the main check:
			)
			{
				++first_shift;
				first_ignored_positions.erase( first_ignored_positions.begin());  // erase first of the remaining erased positions...and go nuts
				++local_first_shift;
			}
		}
		if( msa_b_id < max)
		{
			while( second_ignored_positions.size() > 0 && second_ignored_positions[0] <= msa_b_id)
			{
				++second_shift;
				second_ignored_positions.erase( second_ignored_positions.begin());
				++local_second_shift;
			}
		}

		if( msa_a_id == max && msa_b_id == max)
		{
//			std::cout << "2 gaps in AAMSA by alignment id" << "\n";
			continue;
		}


		// TODO: should one align residues where both AMSAs contain erased residues?
		// add missing information to alignment => if one AMSA has an erased position, assign a gap to the other
		for( int i = 0; i < local_first_shift; ++i)
		{
			tmp = msa_a_id + i + previous_first_shift;
			if( first_sequence[ tmp].GetType() != 'X' && first_sequence[ tmp].GetType() != '-')
			{
				aligned_list.push_back( std::make_pair( tmp, max-1));
			}
		}
		for( int i = 0; i < local_second_shift; ++i)
		{
			tmp = msa_b_id + i + previous_second_shift;
			if( second_sequence[ tmp].GetType() != 'X' && second_sequence[ tmp].GetType() != '-')
			{
				aligned_list.push_back( std::make_pair( max-1, tmp));
			}
		}

		if( msa_a_id != max)
		{
//			std::cout << "shift first msa id from " << msa_a_id << " to " << msa_a_id + first_shift << "\n";
			msa_a_id += first_shift;
			previous_first_shift = first_shift;
		}
		if( msa_b_id != max)
		{
//			std::cout << "shift second msa id from " << msa_b_id << " to " << msa_b_id + second_shift << "\n";
			msa_b_id += second_shift;
			previous_second_shift = second_shift;
		}

		if( msa_a_id == max && ( second_sequence[ msa_b_id].GetType() == 'X' || second_sequence[ msa_b_id].GetType() == '-'))
		{
//			std::cout << "2 gaps by alignment id and amino acid type" << "\n";
			continue;
		}
		if( msa_b_id == max && ( first_sequence[ msa_a_id].GetType() == 'X' || first_sequence[ msa_a_id].GetType() == '-'))
		{
//			std::cout << "2 gaps by amino acid type and alignment id" << "\n";
			continue;
		}
		if( msa_a_id != max && msa_b_id != max
				&& ( first_sequence[ msa_a_id].GetType() == 'X' || first_sequence[ msa_a_id].GetType() == '-')
				&& ( second_sequence[ msa_b_id].GetType() == 'X' || second_sequence[ msa_b_id].GetType() == '-'))
		{
//			std::cout << "2 gaps by amino acid type and alignment id" << "\n";
			continue;
		}

#ifdef DEBUG
		if( msa_a_id < max && msa_b_id < max)
		{
			std::cout << "push back: " << msa_a_id << "   " << first_sequence[ msa_a_id].GetType() << "  " <<  msa_b_id << "  " << second_sequence[ msa_b_id].GetType() << "\n";
		}
		else if( msa_a_id < max)
		{
			std::cout << "push back: " << msa_a_id << "   " << first_sequence[ msa_a_id].GetType() << "  " <<  msa_b_id << "  GAP" << "\n";
		}
		else
		{
			std::cout << "push back: " << msa_a_id << "  GAP  " << "  " <<  msa_b_id << "  " << second_sequence[ msa_b_id].GetType() << "\n";
		}
#endif

		aligned_list.push_back( std::make_pair( msa_a_id, msa_b_id));

	}


	// add last pieces of the sequences that were not part of the AMSA (averaged MSA)
	int count = first_sequence.size() - 1;
	std::vector< std::pair< int,int> > local_rubber;
	for( unsigned int i = 0; i < first_ignored_positions.size(); ++i, --count)
	{
		if( first_sequence[ count].GetType() != '-' && first_sequence[ count].GetType() != 'X')
		{
			local_rubber.insert( local_rubber.begin(), std::make_pair( count, max));
		}
	}
	for( std::vector< std::pair< int,int> >::const_iterator ctr = local_rubber.begin(); ctr !=  local_rubber.end(); ++ctr)
	{
		aligned_list.push_back( *ctr);
	}

	count = second_sequence.size() - 1;
	local_rubber.clear();
	for( unsigned int i = 0; i < second_ignored_positions.size(); ++i, --count)
	{
		if( second_sequence[ count].GetType() != '-' && second_sequence[ count].GetType() != 'X')
		{
			local_rubber.insert( local_rubber.begin(), std::make_pair( max, count));
		}
	}
	for( std::vector< std::pair< int,int> >::const_iterator ctr = local_rubber.begin(); ctr !=  local_rubber.end(); ++ctr)
	{
		aligned_list.push_back( *ctr);
	}



	return std::make_pair( ALIGNED_AMSA.first, aligned_list);
}



inline
std::pair< double, std::vector< std::pair< int, int> > >
ExtractPairAlignmentFromAlignedMSA
(
	const std::pair< double, std::vector< std::pair< int, int> > >  &ALIGNED_AMSA,
	const std::vector< Sequence> &FIRST_MSA,
	const int FIRST_ID,
	const std::vector< int> &FIRST_IGNORED_POSITIONS,
	const std::vector< Sequence> &SECOND_MSA,
	const int SECOND_ID,
	const std::vector< int> &SECOND_IGNORED_POSITIONS
)
{
	return ExtractPairAlignmentFromAlignedMSA( ALIGNED_AMSA, FIRST_MSA[ FIRST_ID], FIRST_IGNORED_POSITIONS, SECOND_MSA[ SECOND_ID], SECOND_IGNORED_POSITIONS);
}




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
