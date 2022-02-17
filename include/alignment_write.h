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
//!  alignment_write.h
//!
//! A collection of functions for writing alignments.
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
//! @author: Lucy Forrest update version number 2021.03.05
/////////////////////////////////////////////////////////////////////////


#ifndef ALIGNMENT_WRITE_H
#define ALIGNMENT_WRITE_H

#include <cstdlib>
#include "command_line_manager.h"
#include "triplet.t.h"
#include "sequence.h"
#include "alignment_variables.h"


inline
std::pair< double, std::vector< std::pair< int, int> > >
ExtractPairAlignmentFromAlignedMSA
(
	const std::pair< double, std::vector< std::pair< int, int> > >  &ALIGNED_AMSA,
	const Sequence &FIRST_MSA_SEQ,
	const std::vector< int> &FIRST_IGNORED_POSITIONS,
	const Sequence &SECOND_MSA_SEQ,
	const std::vector< int> &SECOND_IGNORED_POSITIONS
);


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
);





void WriteAlignmentIndices( const std::pair< double, std::vector< std::pair< int, int> > > &ALIGNMENT, const std::string &FILE, const size_t &LINE_LENGTH);



void WriteAlignedSequencesInClustalwFormat
(
		const std::pair< double, std::vector< std::pair< int, int> > > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND,
		const std::string &FIRST_FASTA_HEADER,
		const std::string &SECOND_FASTA_HEADER,
		const std::string &FILE,
		const size_t &LINE_LENGTH,
		const double &GAP_EXTENSION_PENALTY,
		const std::vector< Triplet<int,int,double> > &ANCHORS,
		const bool WRITE_CLUSTAL_HEADER = true
);

void WriteAlignedSequencesInFastaFormat
(
		const std::pair< double, std::vector< std::pair< int, int> > > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND,
		const std::string &FIRST_FASTA_HEADER,
		const std::string &SECOND_FASTA_HEADER,
		const std::string &FILE,
		const size_t &LINE_LENGTH,
		const double &GAP_EXTENSION_PENALTY
);

void WriteAlignedProfiles
(
		const std::pair< double, std::vector< std::pair< int, int> > > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND,
		const std::string &FILE,
		const double &GAP_EXTENSION_PENALTY,
		const std::string &GAP_VALUE,
		const double &LAST_ELEMENT,
		const ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORES,
		const std::vector< Triplet<int,int,double> > &ANCHORS,
		const std::vector<std::string> &WRITE_PROFILES_HEADER1 = std::vector<std::string>(),
		const std::vector<std::string> &WRITE_PROFILES_HEADER2 = std::vector<std::string>()
);





std::vector< std::pair< int, int> >
AlignmentToSequenceIDs
(
		const std::vector< std::pair< int, int> > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND
);



void WriteAlignment( CommandLineManager &CMD, AlignmentVariables &VARS);


#endif // ALIGNMENT_WRITE_H
