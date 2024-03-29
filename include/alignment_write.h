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
#include "anchor.h"

void WriteAlignmentIndices( const std::pair< double, std::vector< std::pair< int, int> > > &ALIGNMENT, const std::string &FILE, const size_t &LINE_LENGTH)
{
  std::ofstream write( FILE.c_str());

  size_t
    nr_blocks( ALIGNMENT.second.size() / LINE_LENGTH);  // devision of int's !!

  size_t
    remaining( nr_blocks % LINE_LENGTH);

  std::vector< std::pair< int, int> > alignment;
  std::copy( ALIGNMENT.second.begin(), ALIGNMENT.second.end(), alignment.begin());
  std::vector< std::pair< int, int> >::const_iterator itr;

  write << "score: " <<  ALIGNMENT.first << "\n";



  for( size_t b = 0; b < nr_blocks; ++b)
    {

      itr = alignment.begin();
      itr = itr + int( b * LINE_LENGTH);
      for( size_t i = 0; i < LINE_LENGTH; ++i, ++itr)
	{
	  if( itr->first == std::numeric_limits< int>::max())
	    {
	      write << "- ";
	    }
	  else
	    {
	      write << itr->first << " ";
	    }
	}
      write << "\n";

      itr = alignment.begin();
      itr += b * LINE_LENGTH;
      for( size_t i = 0; i < LINE_LENGTH; ++i, ++itr)
	{
	  if( itr->second == std::numeric_limits< int>::max())
	    {
	      write << "- ";
	    }
	  else
	    {
	      write << itr->second << " ";
	    }
	}
      write << "\n";
      write << "\n";
      write << "\n";
    }

  
  for( size_t i = 0; i < remaining; ++i, ++itr)
    {
      if( itr->first == std::numeric_limits< int>::max())
	{
	  write << "- ";
	}
      else
	{
	  write << itr->first << " ";
	}
    }
  write << "\n";
  itr = alignment.begin() + nr_blocks * LINE_LENGTH;
  for( size_t i = 0; i < remaining; ++i, ++itr)
    {

	  if( itr->second == std::numeric_limits< int>::max())
	    {
	      write << "- ";
	    }
	  else
	    {
	      write << itr->second << " ";
	    }

    }

  write.close();
  write.clear();
}



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
)
{
  std::ofstream
	  write( FILE.c_str());

  double
	  identity = 0;

  int
	  non_gaps_in_first_seq = 0,
	  non_gaps_in_second_seq = 0;

  size_t
    nr_blocks( ALIGNMENT.second.size() / LINE_LENGTH);  // devision of int's !!

  size_t
    remaining( ALIGNMENT.second.size() % LINE_LENGTH);

  std::string
	  first_header = FIRST_FASTA_HEADER,
	  second_header = SECOND_FASTA_HEADER;


  if( WRITE_CLUSTAL_HEADER)
  {
	  write << "CLUSTAL W formatted alignment obtained with AlignMe 1.2.2\n\n";
  }

  first_header.resize( 13, ' ');
  second_header.resize( 13, ' ');

  for( size_t i = 0; i < 13; ++i)
  {
	  if (first_header.substr(i,1) == " ")
	  {
		  first_header.replace(i,1, "_");
	  }
	  if (second_header.substr(i,1) == " ")
	  {
		  second_header.replace(i,1, "_");
	  }
  }

  first_header += "   ";
  second_header += "   ";

  DebugWrite ( "nr-blocks: " << nr_blocks);
  DebugWrite ( "remaining: " << remaining);

  std::vector< std::pair< int, int> > alignment( ALIGNMENT.second.size());
  std::copy( ALIGNMENT.second.begin(), ALIGNMENT.second.end(), alignment.begin());
  std::vector< std::pair< int, int> >::const_iterator itr;


  for( size_t b = 0; b < nr_blocks; ++b)
  {

      itr = alignment.begin() + b * LINE_LENGTH;
      write << first_header;
      for( size_t i = 0; i < LINE_LENGTH; ++i, ++itr)
      {
    	  if( itr->first == std::numeric_limits< int>::max())
    	  {
    		  write << "-";
    	  }
    	  else if( itr->first == std::numeric_limits< int>::max()-1)
    	  {
    		  write << ".";
    	  }
    	  else
    	  {
		      write << FIRST[ itr->first].GetType();
    	  }
      }
      write << "\n";

      itr = alignment.begin() + b * LINE_LENGTH;
      write << second_header;
      for( size_t i = 0; i < LINE_LENGTH; ++i, ++itr)
      {
		  if( itr->second == std::numeric_limits< int>::max())
		  {
			  write << "-";
		  }
		  else if( itr->second == std::numeric_limits< int>::max() - 1)
		  {
			  write << ".";
		  }
		  else
		  {
			  write << SECOND[ itr->second].GetType();
		  }
      }
      write << "\n";
      write << "                ";
      itr = alignment.begin() + b * LINE_LENGTH;
      for( size_t i = 0; i < LINE_LENGTH; ++i, ++itr)
      {
    	  if (   itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max() && itr->first != std::numeric_limits< int>::max() -1 && itr->second != std::numeric_limits< int>::max() -1 )
    	  {
    		  if( Contains( ANCHORS, itr->first+1, itr->second+1) )
    		  {
    			  write << "a";
    		  }
    		  else if( FIRST[ itr->first].GetType() == SECOND[ itr->second].GetType())
    		  {
    			  write << "*";
    		  }
    		  else
    		  {
    			  write << " ";
    		  }
    	 }
    	 else
    	 {
				 write << " ";
    	 }
      }

      write << "\n";
      write << "\n";
      write << "\n";
  }


  itr = alignment.begin() + nr_blocks * LINE_LENGTH;
  write << first_header;
  for( size_t i = 0; i < remaining; ++i, ++itr)
  {
      if( itr->first == std::numeric_limits< int>::max())
      {
    	  write << "-";
      }
	  else if( itr->first == std::numeric_limits< int>::max()-1)
	  {
		  write << ".";
	  }
      else
      {
    	  write <<FIRST[  itr->first].GetType();
      }
  }

  write << "\n";
  itr = alignment.begin() + nr_blocks * LINE_LENGTH;
  write << second_header;
  for( size_t i = 0; i < remaining; ++i, ++itr)
    {

	  if( itr->second == std::numeric_limits< int>::max())
	    {
	      write << "-";
	    }
	  else if( itr->second == std::numeric_limits< int>::max()-1)
	  {
		  write << ".";
	  }
	  else
	    {
	      write << SECOND[ itr->second].GetType();
	    }

    }
  write << "\n";

  write << "                ";

  itr = alignment.begin() + nr_blocks * LINE_LENGTH;
  for( size_t i = 0; i < remaining; ++i, ++itr)
  {
	  if (   itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max() && itr->first != std::numeric_limits< int>::max() -1 && itr->second != std::numeric_limits< int>::max() -1 )
	  {
		  if( Contains( ANCHORS, itr->first+1, itr->second+1) )
		  {
			  write << "a";
		  }
		  else if( FIRST[ itr->first].GetType() == SECOND[ itr->second].GetType())
		  {
			  write << "*";
		  }
		  else
		  {
			  write << " ";
		  }
	 }
	 else
	 {
			 write << " ";
	 }
  }
  write << "\n";
  write << "\n";
  for( itr = alignment.begin(); itr != alignment.end(); ++itr)
  {
	  if( itr->first != std::numeric_limits< int>::max() ||  itr->first != std::numeric_limits< int>::max() - 1)
	  {
		  non_gaps_in_first_seq++;
	  }
	  if (itr->second != std::numeric_limits< int>::max() ||  itr->second != std::numeric_limits< int>::max() - 1)
	  {
		  non_gaps_in_second_seq++;
	  }
  }


  for( itr = alignment.begin(); itr != alignment.end(); ++itr)
  {

	  if
	  (
			  itr->first != std::numeric_limits< int>::max()
			  && itr->second != std::numeric_limits< int>::max()
			  && itr->first != std::numeric_limits< int>::max() - 1
			  && itr->second != std::numeric_limits< int>::max() - 1
			  && FIRST[ itr->first].GetType() == SECOND[ itr->second].GetType()
	  )
	  {
		  identity += 1.0;
	  }
  }

  if( non_gaps_in_first_seq > non_gaps_in_second_seq)
  {
	  identity /= double (non_gaps_in_second_seq);
  }
  else
  {
	  identity /= double (non_gaps_in_first_seq);
  }
  identity *= 100;

  write.close();
  write.clear();
}


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
)
{
  std::ofstream
	  write( FILE.c_str());

  double
	  identity = 0;

  int
	  non_gaps_in_first_seq = 0,
	  non_gaps_in_second_seq = 0;

  size_t
    nr_blocks( ALIGNMENT.second.size() / LINE_LENGTH);  // devision of int's !!

  size_t
    remaining( ALIGNMENT.second.size() % LINE_LENGTH);

  std::string
	  first_header = FIRST_FASTA_HEADER,
	  second_header = SECOND_FASTA_HEADER;

  first_header.resize( 13, ' ');
  second_header.resize( 13, ' ');

  first_header += "";
  second_header += "";

  DebugWrite ( "nr-blocks: " << nr_blocks);
  DebugWrite ( "remaining: " << remaining);

  std::vector< std::pair< int, int> > alignment( ALIGNMENT.second.size());
  std::copy( ALIGNMENT.second.begin(), ALIGNMENT.second.end(), alignment.begin());
  std::vector< std::pair< int, int> >::const_iterator itr;

  write << ">" << first_header << "\n";
  for( size_t b = 0; b < nr_blocks; ++b)
  {

      itr = alignment.begin() + b * LINE_LENGTH;
      for( size_t i = 0; i < LINE_LENGTH; ++i, ++itr)
      {
    	  if( itr->first == std::numeric_limits< int>::max())
    	  {
    		  write << "-";
    	  }
    	  else if( itr->first == std::numeric_limits< int>::max()-1)
    	  {
    		  write << ".";
    	  }
    	  else
    	  {
		      write << FIRST[ itr->first].GetType();
    	  }
      }
      write << "\n";
  }
  itr = alignment.begin() + nr_blocks * LINE_LENGTH;
  for( size_t i = 0; i < remaining; ++i, ++itr)
  {
      if( itr->first == std::numeric_limits< int>::max())
      {
    	  write << "-";
      }
	  else if( itr->first == std::numeric_limits< int>::max()-1)
	  {
		  write << ".";
	  }
      else
      {
    	  write <<FIRST[  itr->first].GetType();
      }
  }
  write << "\n";


  write << ">" << second_header << "\n";
  for( size_t b = 0; b < nr_blocks; ++b)
  {
      itr = alignment.begin() + b * LINE_LENGTH;
      for( size_t i = 0; i < LINE_LENGTH; ++i, ++itr)
      {
		  if( itr->second == std::numeric_limits< int>::max())
		  {
			  write << "-";
		  }
		  else if( itr->second == std::numeric_limits< int>::max() - 1)
		  {
			  write << ".";
		  }
		  else
		  {
			  write << SECOND[ itr->second].GetType();
		  }
      }
      write << "\n";
  }

  itr = alignment.begin() + nr_blocks * LINE_LENGTH;
  for( size_t i = 0; i < remaining; ++i, ++itr)
    {

	  if( itr->second == std::numeric_limits< int>::max())
	    {
	      write << "-";
	    }
	  else if( itr->second == std::numeric_limits< int>::max()-1)
	  {
		  write << ".";
	  }
	  else
	    {
	      write << SECOND[ itr->second].GetType();
	    }

    }
  write << "\n";

  for( itr = alignment.begin(); itr != alignment.end(); ++itr)
  {
	  if( itr->first != std::numeric_limits< int>::max() ||  itr->first != std::numeric_limits< int>::max() - 1)
	  {
		  non_gaps_in_first_seq++;
	  }
	  if (itr->second != std::numeric_limits< int>::max() ||  itr->second != std::numeric_limits< int>::max() - 1)
	  {
		  non_gaps_in_second_seq++;
	  }
  }


  for( itr = alignment.begin(); itr != alignment.end(); ++itr)
  {

	  if
	  (
			  itr->first != std::numeric_limits< int>::max()
			  && itr->second != std::numeric_limits< int>::max()
			  && itr->first != std::numeric_limits< int>::max() - 1
			  && itr->second != std::numeric_limits< int>::max() - 1
			  && FIRST[ itr->first].GetType() == SECOND[ itr->second].GetType()
	  )
	  {
		  identity += 1.0;
	  }
  }

  if( non_gaps_in_first_seq > non_gaps_in_second_seq)
  {
	  identity /= double (non_gaps_in_second_seq);
  }
  else
  {
	  identity /= double (non_gaps_in_first_seq);
  }
  identity *= 100;

  write.close();
  write.clear();
}


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
)
{
//	std::cout << __FUNCTION__ << " " << GAP_EXTENSION_PENALTY << "  " << FILE << std::endl;
  std::ofstream write( FILE.c_str());
  if( !write)
  {
	  std::cout << "WARNING: It was not possible to open the file (" << FILE << ") to which the profiles should be written to! \n";
  }

  std::vector< std::pair< int, int> > alignment( ALIGNMENT.second.size());
  std::copy( ALIGNMENT.second.begin(), ALIGNMENT.second.end(), alignment.begin());
  std::vector< std::pair< int, int> >::const_iterator itr = alignment.begin();

  std::vector< double>
	  first_vector,
	  second_vector;

// bool
//	  write_prof_dep_sub_matrix = false,
//	  write_sub_matrix = false,
//	  write_pssm = false;

  size_t
	  number_profiles( FIRST[0].GetProfiles().size()),
	  cc( 0);

  DebugWrite("First profiles size " << FIRST[0].GetProfiles().size());
  DebugWrite("Second profiles size " << SECOND[0].GetProfiles().size());

  cc = 1;
  write << "#Column  "  << ":  the position within the alignment\n";

  SumFunction< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> *
	  sum = ( SumFunction< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> *) SCORES.get();
  std::vector< std::pair< double, ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > > >
	  functions = sum->GetData();
  std::vector< std::pair< double, ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > > >::const_iterator
	  atr;
  ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> >
	  prof_sub, sub, pssm;

  // WRITE PROCESS FOR MATRICES IS COMMENTED OUT FOR THE WEBSERVER!
  /*
  for( atr = functions.begin(); atr != functions.end(); ++atr)
  {
	  if( atr->second->GetClassID() == ScoreProfileDependentSequenceSimilarity().GetClassID())
	  {
		  write << "#Column " << cc++ << " is the substitution matrix profile " << std::endl;
		  prof_sub = atr->second;
		  write_prof_dep_sub_matrix = true;
	  }
	  if( atr->second->GetClassID() == ScoreSequenceSimilarity().GetClassID())
	  {
		  write << "#Column " << cc++ << " is the substitution matrix profile " << std::endl;
		  sub = atr->second;
		  write_sub_matrix = true;
	  }
	  if( atr->second->GetClassID() == ScorePositionSpecificSimilarity().GetClassID())
	  {
		  write << "#Column " << cc++ << " is the position specific substitution matrix profile" << std::endl;
		  pssm = atr->second;
		  write_pssm = true;
	  }
  }
*/
  for (size_t i = 0; i < WRITE_PROFILES_HEADER1.size(); ++i)
  {
	  write << "#Column " << std::right << std::setw(2) << c++ << ":  " << WRITE_PROFILES_HEADER1[i];
  }

  for (size_t i = 0; i < WRITE_PROFILES_HEADER2.size(); ++i)
  {
	  write << "#Column " << std::right << std::setw(2) << c++ << ":  " << WRITE_PROFILES_HEADER2[i];
  }
  write << "#Column " << std::right << std::setw(2) << c++ << ":  residue id within the first sequence\n";
  write << "#Column " << std::right << std::setw(2) << c++ << ":  one letter amino acid type of the first sequence\n";
  write << "#Column " << std::right << std::setw(2) << c++ << ":  residue id within the second sequence\n";
  write << "#Column " << std::right << std::setw(2) << c++ << ":  one letter amino acid type of the second sequence\n";
  write << "#Column " << std::right << std::setw(2) << c++ << ":  conserved positions ('*')\n";
  write << "#Column " << std::right << std::setw(2) << c++ << ":  anchors for the first sequence  ('a')\n";
  write << "#Column " << std::right << std::setw(2) << c++ << ":  anchors for the second sequence ('a')\n";
  write << "#Conservation and anchors are shown only if applicable" << std::endl;
  write << "#Gaps are represented by: "<< GAP_VALUE << "\n  " << std::endl;

  cc = 0;
  for( ; itr != alignment.end(); ++itr, ++cc)
  {
	  write.width(4);
	  write << cc << "  ";
/*
	  if( itr->first != std::numeric_limits< int>::max() &&  itr->second != std::numeric_limits< int>::max())
	  {
		  if( write_prof_dep_sub_matrix)
		  {
			  write.width( 14);
			  write << (*prof_sub)( std::make_pair( FIRST[ itr->first], SECOND[ itr->second])) << "  ";
//			  write.flush();
		  }
		  if( write_sub_matrix)
		  {
			  write.width( 14);
			  write << (*sub)( std::make_pair( FIRST[ itr->first], SECOND[ itr->second])) << "  ";
//			  write.flush();
		  }
		  if( write_pssm)
		  {
			  write.width( 14);
			  write << (*pssm)( std::make_pair( FIRST[ itr->first], SECOND[ itr->second])) << "  ";
//			  write.flush();
		  }
	  }
	  else
	  {
		  if( write_prof_dep_sub_matrix)
		  {
			  write.width( 14);
			  write << GAP_VALUE << "  ";
//			  write.flush();
		  }
		  if( write_sub_matrix)
		  {
			  write.width( 14);
			  write << GAP_VALUE << "  ";
//			  write.flush();
		  }
		  if( write_pssm)
		  {
			  write.width( 14);
			  write << GAP_VALUE << "  ";
//			  write.flush();
		  }
	  }
*/
	  if( itr->first != std::numeric_limits< int>::max())
	  {
		  first_vector = FIRST[ itr->first].GetProfiles();
		  DebugWriteNoFlush( FIRST[ itr->first].GetType());
		  for( std::vector< double>::const_iterator vitr = first_vector.begin(); vitr != first_vector.end(); ++vitr)
		  {
			  write.setf( std::ios_base::fixed, std::ios_base::floatfield );
			  write.precision(3);
			  write.width( 8);
			  write << *vitr << "  ";
		  }
	  }
	  else
	  {
		  DebugWriteNoFlush( "gap");
		  for( size_t i = 0; i < number_profiles; ++i)
		  {
			  write.width( 8);
			  write << GAP_VALUE << "  ";
		  }
	  }
	  if( itr->second != std::numeric_limits< int>::max())
	  {
		  second_vector = SECOND[ itr->second].GetProfiles();
		  DebugWrite( SECOND[ itr->second].GetType());
		  for( std::vector< double>::const_iterator vitr = second_vector.begin(); vitr != second_vector.end(); ++vitr)
		  {
			  write.precision(3);
			  write.width( 8);
			  write << *vitr << "  ";
		  }
	  }
	  else
	  {
		  DebugWriteNoFlush( "gap");
		  for( size_t i = 0; i < number_profiles; ++i)
		  {
			  write.width( 8);
			  write << GAP_VALUE << "  ";
		  }
	  }

	  if( itr->first != std::numeric_limits< int>::max())
	  {
		  write.width(6);
		  write << itr->first+1 << " ";
		  write.width( 3);
		  write << FIRST[ itr->first].GetType() << " ";
	  }
	  else
	  {
		  write.width( 6);
		  write << GAP_VALUE << " ";
		  write.width( 3);
		  write << GAP_VALUE << " ";
	  }
	  if( itr->second != std::numeric_limits< int>::max())
	  {
		  write.width( 6);
		  write << itr->second+1 << " ";
		  write.width( 3);
		  write << SECOND[ itr->second].GetType() << " ";
	  }
	  else
	  {
		  write.width( 6);
		  write << GAP_VALUE << " ";
		  write.width( 3);
		  write << GAP_VALUE << " ";
	  }
	  if(  itr->first != std::numeric_limits< int>::max() && itr->second != std::numeric_limits< int>::max() && FIRST[ itr->first].GetType() == SECOND[ itr->second].GetType() )
	  {
		  write << " *";
	  }
	  else
	  {
		  write << "  ";
	  }
	  if( ContainsFirst( ANCHORS, itr->first+1 ) )
	  {
		  write << "  a   ";
	  }
	  else
	  {
		  write << "      ";
	  }
	  if( ContainsSecond( ANCHORS, itr->second+1 ) )
	  {
		  write << "a";
	  }

	  write << "" << std::endl;
  }

  write << "" << std::endl;
  // VERSION 1.1
//  write.setf( std::ios::showpoint | std::ios::fixed, std::ios::floatfield);
//  write.precision( 6);
//  write << "#normalized_DM_last_element: " << LAST_ELEMENT / NrAligned( ALIGNMENT.second) << "\n";
//  write << "#normalized_last_element: " << LAST_ELEMENT / LengthAligned( ALIGNMENT.second) << "\n";
  write.close();
  write.clear();
}







std::vector< std::pair< int, int> >
AlignmentToSequenceIDs
(
		const std::vector< std::pair< int, int> > &ALIGNMENT,
		const Sequence &FIRST,
		const Sequence &SECOND
)
{
	std::vector< std::pair< int, int> >
		list;

	int
		first_count = 0,
		second_count = 0;



	for( std::vector< std::pair< int, int> >::const_iterator itr = ALIGNMENT.begin(); itr != ALIGNMENT.end(); ++itr)
	{
		if( itr->first  < std::numeric_limits< int>::max() - 10  &&
    		  itr->second < std::numeric_limits< int>::max() - 10  &&
    		  SECOND[ itr->second].GetType() != '-' &&  SECOND[ itr->second].GetType() != '=' &&
    		  FIRST[ itr->first].GetType() != '-' &&  FIRST[ itr->first].GetType() != '=')
			{
				list.push_back( std::make_pair( first_count, second_count));
				++first_count;
				++second_count;
//				std::cout << "pb: " << first_count << "  " << second_count << " type: " << FIRST[ itr->first].GetType() << " " << SECOND[ itr->second].GetType() <<  std::endl;
			}
		else
		{
			if (itr->first  != std::numeric_limits< int>::max() && itr->first  != std::numeric_limits< int>::max()-1)
			{
				if( FIRST[ itr->first].GetType() != '-' &&  FIRST[ itr->first].GetType() != '=')
				{
					++first_count;
				}
			}

			if (itr->second  != std::numeric_limits< int>::max() && itr->second  != std::numeric_limits< int>::max()-1)
			{
				if( SECOND[ itr->second].GetType() != '-' &&  SECOND[ itr->second].GetType() != '=')
				{
					++second_count;
				}
			}
		}
	}
  return list;
}







void WriteAlignment( CommandLineManager &CMD, AlignmentVariables &VARS)
{
	  if (VARS.outputfile_aligned_sequences != "" && VARS.first_fasta_file != "" && VARS.second_fasta_file != "")
	  {
		  if ( VARS.alignment_output_format == "clustalw")
		  {
			  DebugWrite( "write aligned sequences in clustalw format");
			  WriteAlignedSequencesInClustalwFormat( VARS.score_and_alignment, VARS.first_sequence, VARS.second_sequence, VARS.first_fasta_id, VARS.second_fasta_id, VARS.outputfile_aligned_sequences, VARS.line_length, VARS.gap_extension_penalty, VARS.anchors);
			  DebugWrite( "done");
		  }
		  else if ( VARS.alignment_output_format == "fasta")
		  {
			  DebugWrite( "write aligned sequences in fasta format");
			  WriteAlignedSequencesInFastaFormat( VARS.score_and_alignment, VARS.first_sequence, VARS.second_sequence, VARS.first_fasta_id, VARS.second_fasta_id, VARS.outputfile_aligned_sequences, VARS.line_length, VARS.gap_extension_penalty);
			  DebugWrite( "done");
		  }
		  else
		  {
			  DebugWrite( "write aligned sequences in clustalw format");
			  WriteAlignedSequencesInClustalwFormat( VARS.score_and_alignment, VARS.first_sequence, VARS.second_sequence, VARS.first_fasta_id, VARS.second_fasta_id, VARS.outputfile_aligned_sequences, VARS.line_length, VARS.gap_extension_penalty, VARS.anchors);
			  DebugWrite( "done");
		  }
	  }
	  else if (VARS.outputfile_aligned_sequences != "")
	  {
		  std::cout << "WARNING: The flag -outputfile_aligned_sequences is provided but not the flags -fasta_file1 and -fasta_file2. However, a pairwise sequence alignment can only be generated if two sequences in fasta format are provided. Because this is not the case, such an output is not created! \n";
	  }

	  if( VARS.outputfile_aligned_profiles != "")
	  {
		  DebugWrite( "write aligned profiles");
		  if( CMD.IsFlagSet( "profile_gap_value_for_plotting"))
		  {
			  VARS.gap_value =  CMD.GetFirstArgument( "profile_gap_value_for_plotting");
		  }
		  WriteAlignedProfiles( VARS.score_and_alignment, VARS.first_sequence, VARS.second_sequence,  VARS.outputfile_aligned_profiles , VARS.gap_extension_penalty, VARS.gap_value, VARS.last_element, VARS.scores, VARS.anchors, VARS.write_profiles_header1, VARS.write_profiles_header2 );
		  DebugWrite( "done");
	  }

	  if( CMD.IsFlagSet( "output_extracted_sequences") && !CMD.IsFlagSet( "extract_from_MSA_sequences_with_ids"))
	  {
		  std::cout << "WARNING: To use the flag -output_extracted_sequences, the flag -extract_from_MSA_sequences_with_ids has to be set. This is not the case and therefore, sequences are not extracted from the averaged MSAs! \n";

	  }
	  if( CMD.IsFlagSet( "extract_from_MSA_sequences_with_ids"))
	  {
		  if (VARS.first_msa_file != "" && VARS.second_msa_file != "")
		  {
			  std::vector< size_t>
				  ids = CMD.GetConvertedArguments< size_t>( "extract_from_MSA_sequences_with_ids");

#ifdef SECURE

			  if (ids.size() != 2)
			  {
				  std::cerr << "ERROR: After the flag -extract_from_MSA_sequences_with_ids, two numbers have to be provided! \n";
				  exit(-1);
			  }

				if (VARS.first_msa.size() <= ids[0])
				{
					std::cout << "ERROR: From msa_file1 sequence #"<<  ids[0] <<" can not be extracted because there are only "<< VARS.first_msa.size()  <<" sequences in msa_file1!\n";
					exit(-1);
				}
				if ( VARS.second_msa.size() <= ids[1])
				{
					std::cout << "ERROR: From msa_file2 sequence #"<<  ids[1] <<" can not be extracted because there are only "<< VARS.first_msa.size()  <<" sequences in msa_file2!\n";
					exit(-1);
				}
#endif

			  std::pair< double, std::vector< std::pair< int, int> > >
				aligned_pair = ExtractPairAlignmentFromAlignedMSA( VARS.score_and_alignment, VARS.first_msa, ids[0], VARS.first_ignored_positions, VARS.second_msa, ids[1], VARS.second_ignored_positions);

		//	  std::ofstream write;
		//	  write.open( "aligned_pairs.txt");
		//	  write << aligned_pair.second;
		//	  write.close();

			 // DebugWrite( "write aligned sequences");

			  if( CMD.IsFlagSet( "output_extracted_sequences"))
			  {
				  VARS.outputfile_aligned_sequences = CMD.GetFirstArgument( "output_extracted_sequences");
			  }
			  else
			  {
				  VARS.outputfile_aligned_sequences = "output_extracted_sequences_from_MSA.txt";
				  std::cout << "WARNING: The flag -output_extracted_sequences is missing although the flag -extract_from_MSA_sequences_with_ids is provided. The extracted sequences will be written to " << VARS.outputfile_aligned_sequences << std::endl;
			  }

			  // CHECK OUTPUTFILE_ALIGNED_IDS


			  if ( VARS.alignment_output_format == "clustalw")
			  {
				  WriteAlignedSequencesInClustalwFormat( aligned_pair, VARS.first_msa[ ids[ 0]], VARS.second_msa[ ids[ 1]], "seq_of_msa1", "seq_of_msa2", VARS.outputfile_aligned_sequences, VARS.line_length, VARS.gap_extension_penalty, VARS.anchors, false); //60 is just for printing the alignment
			  }
			  else if ( VARS.alignment_output_format == "fasta")
			  {
				  WriteAlignedSequencesInFastaFormat( aligned_pair, VARS.first_msa[ ids[ 0]], VARS.second_msa[ ids[ 1]], "seq_of_msa1", "seq_of_msa2", VARS.outputfile_aligned_sequences, VARS.line_length, VARS.gap_extension_penalty); //60 is just for printing the alignment
			  }
			  else
			  {
				  WriteAlignedSequencesInClustalwFormat( aligned_pair, VARS.first_msa[ ids[ 0]], VARS.second_msa[ ids[ 1]], "seq_of_msa1", "seq_of_msa2", VARS.outputfile_aligned_sequences, VARS.line_length, VARS.gap_extension_penalty, VARS.anchors, false); //60 is just for printing the alignment
			  }
		  }
		  else
		  {
			  std::cout << "WARNING: The flag -extract_from_MSA_sequences_with_ids can only be used in combination with the flags -msa_file1 and -msa_file2 which is not the case. Therefore, such an output is not created! \n";
		  }
	  }

}

#endif // ALIGNMENT_WRITE_H
