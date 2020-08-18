/////////////////////////////////////////////////////////////////////////
//
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
//
//!
//!
//! alignment_write.h
//!
//! A collection of functions for writing alignments.
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef ALIGNMENT_WRITE_H
#define ALIGNMENT_WRITE_H

#include <cstdlib>

void WriteAlignmentIndices( const std::pair< double, std::list< std::pair< int, int> > > &ALIGNMENT, const std::string &FILE, const size_t &LINE_LENGTH)
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



void WriteAlignedSequences
(
		const std::pair< double, std::list< std::pair< int, int> > > &ALIGNMENT,
		const AASequence &FIRST,
		const AASequence &SECOND,
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

  first_header += "  ";
  second_header += "  ";

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
    		  write << "=";
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
			  write << "=";
		  }
		  else
		  {
			  write << SECOND[ itr->second].GetType();
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
		  write << "=";
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
		  write << "=";
	  }
	  else
	    {
	      write << SECOND[ itr->second].GetType();
	    }

    }
  write << "\n" << "\n";

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
		const std::pair< double, std::list< std::pair< int, int> > > &ALIGNMENT,
		const AASequence &FIRST,
		const AASequence &SECOND,
		const std::string &FILE,
		const double &GAP_EXTENSION_PENALTY,
		const double &GAP_VALUE,
		const double &LAST_ELEMENT,
		const std::vector<std::string> &WRITE_PROFILES_HEADER1 = std::vector<std::string>(),
		const std::vector<std::string> &WRITE_PROFILES_HEADER2 = std::vector<std::string>()
)
{
  std::ofstream write( FILE.c_str());
  if( !write)
  {
	  std::cout << "WARNING: the output file (" << FILE << ") could not be opened!" << "\n";
  }

  std::vector< std::pair< int, int> > alignment( ALIGNMENT.second.size());
  std::copy( ALIGNMENT.second.begin(), ALIGNMENT.second.end(), alignment.begin());
  std::vector< std::pair< int, int> >::const_iterator itr = alignment.begin();

  std::vector< double>
	  first_vector,
	  second_vector;

  size_t
	  number_profiles( FIRST[0].GetProfiles().size()),
	  cc( 0);

  DebugWrite("First profiles size " << FIRST[0].GetProfiles().size());
  DebugWrite("Second profiles size " << SECOND[0].GetProfiles().size());

  write << "#Column 1 is the position within the alignment \n";

  for (size_t i = 0; i < WRITE_PROFILES_HEADER1.size(); ++i)
  {
	  write << "#Column " << i+2 << WRITE_PROFILES_HEADER1[i];
  }

  for (size_t i = 0; i < WRITE_PROFILES_HEADER2.size(); ++i)
  {
	  write << "#Column " << i+2+WRITE_PROFILES_HEADER1.size() << WRITE_PROFILES_HEADER2[i];
  }

  for( ; itr != alignment.end(); ++itr, ++cc)
  {
	  write << cc << "  ";
	  if( itr->first != std::numeric_limits< int>::max())
	  {
		  first_vector = FIRST[ itr->first].GetProfiles();
		  DebugWriteNoFlush( FIRST[ itr->first].GetType());
		  for( std::vector< double>::const_iterator vitr = first_vector.begin(); vitr != first_vector.end(); ++vitr)
		  {
			  write.width( 14);
			  write << *vitr << "  ";
		  }
	  }
	  else
	  {
		  DebugWriteNoFlush( "gap");
		  for( size_t i = 0; i < number_profiles; ++i)
		  {
			  write.width( 14);
			  write << GAP_VALUE << "  ";
		  }
	  }
	  if( itr->second != std::numeric_limits< int>::max())
	  {
		  second_vector = SECOND[ itr->second].GetProfiles();
		  DebugWrite( SECOND[ itr->second].GetType());
		  for( std::vector< double>::const_iterator vitr = second_vector.begin(); vitr != second_vector.end(); ++vitr)
		  {
			  write.width( 14);
			  write << *vitr << "  ";
		  }
	  }
	  else
	  {
		  DebugWriteNoFlush( "gap");
		  for( size_t i = 0; i < number_profiles; ++i)
		  {
			  write.width( 14);
			  write << GAP_VALUE << "  ";
		  }
	  }
	  write << "\n";
  }

  write << "\n";

  write.close();
  write.clear();
}


#endif // ALIGNMENT_WRITE_H























